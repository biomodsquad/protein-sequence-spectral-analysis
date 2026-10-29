"""Reusable notebook controls for the embedding atlas; no network or fitting at import."""

import hashlib
import json

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform
from sklearn.decomposition import PCA
from sklearn.manifold import trustworthiness
from sklearn.neighbors import NearestNeighbors

from .analyze import COLORS, gram_features, projector_features, shape_distance, unit_profile


def prepare_vectors(catalogue, profiles, families, *, sign_invariant=True,
                    deduplicate=True, decimals=10, axis='both'):
    """Select source families/orientations and build signed or sign-invariant shape features."""
    if axis not in ('both', 'rows', 'columns'):
        raise ValueError("axis must be 'both', 'rows', or 'columns'")
    mask = catalogue.family.isin(families).to_numpy()
    if axis != 'both':
        mask &= catalogue.axis.fillna('').isin(['', axis]).to_numpy()
    rows = catalogue.loc[mask].copy().reset_index(drop=True)
    x = np.asarray(profiles)[mask]
    if len(x) < 3:
        raise ValueError('Choose families/orientations containing at least three channel instances')
    lookup, unique, ids = {}, [], []
    for index, profile in enumerate(x):
        canonical = profile.copy()
        if sign_invariant:
            canonical *= np.sign(canonical[np.argmax(abs(canonical))])
        key = tuple(np.round(canonical, decimals)) if deduplicate else (index,)
        if key not in lookup:
            lookup[key] = len(unique)
            unique.append(profile)
        ids.append(lookup[key])
    if len(unique) < 3:
        raise ValueError('Fewer than three unique shapes remain; expand the source families')
    x = np.asarray(unique)
    rows['unique_id'] = ids
    priority = {'Kidera':0, 'BLOSUM100':1, 'AAindex1':2, 'AAindex2':3, 'AAindex3':4}
    points = (rows.assign(priority=rows.family.map(priority)).sort_values('priority',kind='stable')
              .drop_duplicates('unique_id').set_index('unique_id').loc[np.arange(len(x))]
              .reset_index())
    alias_map = rows.groupby('unique_id').id.apply(list)
    points['alias_count'] = points.unique_id.map(alias_map.map(len))
    points['alias_names'] = points.unique_id.map(alias_map.map(lambda names:', '.join(names)))
    features = projector_features(x) if sign_invariant else x
    similarity, distance = shape_distance(features)
    return rows, points, x, features, similarity, distance


def prepare_whole(representations, families, *, equal_channels=True, axis='both'):
    """Choose whole configurations and compare centered residue Grams with explicit weights."""
    if axis not in ('both','rows','columns'):
        raise ValueError("axis must be 'both', 'rows', or 'columns'")
    chosen = [r for r in representations if r['family'] in families and (
        axis == 'both' or r['orientation'] in (axis,'rows=columns','builtin','all scalars'))]
    if len(chosen) < 3:
        raise ValueError('Choose families/orientations containing at least three whole representations')
    features, records = [], []
    for item in chosen:
        gram, stats = gram_features(item['embedding'], equal_channels)
        features.append(gram)
        records.append({k:v for k,v in item.items() if k != 'embedding'} | stats |
                       {'channels':item['embedding'].shape[1]})
    features = np.array(features)
    similarity, distance = shape_distance(features)
    return pd.DataFrame(records), features, similarity, distance, chosen


def cluster_profiles(distance, thresholds, *, sign_invariant, method='complete'):
    """Cut hierarchical trees in the selected metric; pairwise guarantees require complete linkage."""
    if method not in ('complete','average','single','ward'):
        raise ValueError('Choose complete, average, single, or ward linkage')
    if not thresholds or any(not 0 < value <= 1 for value in thresholds):
        raise ValueError('Similarity thresholds must be in (0,1]')
    tree = linkage(squareform(distance,checks=False),method=method)
    summary, assignments = [], {}
    for value in thresholds:
        cut = np.sqrt(2*(1-(value*value if sign_invariant else value)))
        labels = fcluster(tree,cut,criterion='distance')
        assignments[value] = labels
        summary.append({'cut_similarity':value,'clusters':len(set(labels)),
                        'largest_cluster':int(np.bincount(labels).max()),
                        'pairwise_guarantee':method == 'complete'})
    return tree,pd.DataFrame(summary),assignments


def project_pca(features, components):
    """Return scores and the fitted PCA so users can inspect/load more than two axes."""
    count = min(int(components),len(features)-1,features.shape[1])
    if count < 2:
        raise ValueError('PCA needs at least two components and three observations')
    model = PCA(n_components=count,svd_solver='full')
    return model.fit_transform(features),model


def umap_signature(distance, settings):
    """Fingerprint the actual precomputed distances and effective layout settings."""
    if len(distance) < 4:
        raise ValueError('UMAP needs at least four observations here; include another source family')
    parameters = dict(settings)
    parameters['n_neighbors'] = min(parameters.get('n_neighbors',15),len(distance)-1)
    parameters.update(n_components=2,metric='precomputed',n_jobs=1)
    data = np.asarray(distance,dtype=np.float32)
    signature = hashlib.sha256(data.tobytes()+json.dumps(parameters,sort_keys=True).encode()).hexdigest()
    return data,parameters,signature


def cached_umap(distance, settings, cache, *, use_cache=True):
    """Cache fixed-seed layouts by input distances plus effective UMAP parameters."""
    from umap import UMAP
    data,parameters,signature = umap_signature(distance,settings)
    enabled = use_cache and parameters.get('random_state') is not None
    hit = enabled and signature in cache
    if not hit:
        coordinates = UMAP(**parameters).fit_transform(data)
        if enabled:
            cache[signature] = coordinates
    else:
        coordinates = cache[signature]
    return coordinates,parameters,signature,hit


def scatter(coordinates, table, *, title, x_label='Dimension 1', y_label='Dimension 2',
            marker_size=5, reference_size=11, opacity=.55, height=620, colors=None):
    """Hoverable plot with reference markers, source descriptions, and all available aliases."""
    colors = COLORS if colors is None else colors
    fig = go.Figure()
    for family,color in colors.items():
        indices = np.flatnonzero(table.family.to_numpy() == family)
        if not len(indices):
            continue
        rows = table.iloc[indices]
        reference = family in ('Kidera','BLOSUM100')
        aliases = rows.alias_names if 'alias_names' in rows else rows.id
        hover = [f'{name}<br>{description}<br>Aliases: {other}' for name,description,other in
                 zip(rows.id,rows.description,aliases)]
        fig.add_trace(go.Scattergl(x=coordinates[indices,0],y=coordinates[indices,1],
                                  name=f'{family} (n={len(indices)})',mode='markers',text=hover,
                                  hovertemplate='%{text}<extra></extra>',
                                  marker=dict(color=color,size=reference_size if reference else marker_size,
                                              opacity=1 if reference else opacity,
                                              symbol={'BLOSUM100': 'star', 'Kidera': 'x'}.get(family, 'circle'),
                                              line=dict(color='black', width=0.5 if reference else 0))))
    fig.update_layout(template='plotly_white',height=height,title=title,
                      xaxis_title=x_label,yaxis_title=y_label)
    return fig


def neighbor_table(reference_id, table, similarity, *, count=12):
    """Rank original-space neighbors, independent of the chosen visualization."""
    found = np.flatnonzero(table.id.to_numpy() == reference_id)
    if not len(found):
        raise KeyError(f'{reference_id} not found in the selected source families')
    i = int(found[0])
    indices = np.argsort(-similarity[i],kind='stable')
    indices = indices[indices != i][:count]
    output = table.iloc[indices].copy()
    output.insert(0,'similarity',similarity[i,indices])
    return output


def medoids(table,distance,labels):
    """Choose actual distance medoids; retain groups and original identifiers."""
    representatives = []
    for cluster in sorted(set(labels)):
        indices = np.flatnonzero(labels == cluster)
        choice = indices[np.argmin(distance[np.ix_(indices,indices)].mean(axis=1))]
        row = table.iloc[choice].to_dict()
        row.update(cluster=int(cluster),unique_members=len(indices))
        representatives.append(row)
    return pd.DataFrame(representatives)


def quality(distance, coordinates, *, neighbors=10, audit_size=1000, seed=0):
    """Neighborhood recall and fixed-subset trustworthiness, with explicit self exclusion."""
    n = len(distance)
    audit_size = min(int(audit_size),n)
    k = min(int(neighbors),(audit_size-1)//2)
    if k < 1:
        raise ValueError('Quality audit needs at least three points and a positive neighbor count')
    high = NearestNeighbors(n_neighbors=k+1,metric='precomputed').fit(distance).kneighbors(
        distance,return_distance=False)
    low = NearestNeighbors(n_neighbors=k+1).fit(coordinates).kneighbors(
        coordinates,return_distance=False)
    high = [[j for j in row if j != i][:k] for i,row in enumerate(high)]
    low = [[j for j in row if j != i][:k] for i,row in enumerate(low)]
    recall = float(np.mean([len(set(a)&set(b))/k for a,b in zip(high,low)]))
    rows = np.sort(np.random.default_rng(seed).choice(n,audit_size,replace=False))
    trust = trustworthiness(distance[np.ix_(rows,rows)],coordinates[rows],
                            n_neighbors=k,metric='precomputed')
    return {'neighbors':k,'neighbor_recall':recall,'trustworthiness':float(trust),
            'trustworthiness_audit_size':audit_size,'audit_seed':seed}


def aligned_matrix_features(representations):
    """Signed, label-aligned score features only for matrices with the same canonical axes."""
    chosen = [item for item in representations if item['square']]
    if len(chosen) < 3:
        raise ValueError('At least three square matrix configurations are needed')
    features = np.array([unit_profile(item['embedding'].ravel()) for item in chosen])
    table = pd.DataFrame([{k:v for k,v in item.items() if k != 'embedding'} for item in chosen])
    return table,features


def cka_clustermap(similarity, table, tree, *, colorscale='Viridis',
                   color_range=(0.0, 1.0), height=850, show_labels=False,
                   dendrogram_fraction=0.20):
    """Show original-space CKA with aligned row/column linkage dendrograms.

    ``tree`` must be fitted to sqrt(2 * (1 - CKA)) in the same input order.
    Both axes use its leaf order; branch heights remain linkage distances.
    """
    from scipy.cluster.hierarchy import dendrogram
    from plotly.subplots import make_subplots

    if not 0 < dendrogram_fraction < 0.5:
        raise ValueError('dendrogram_fraction must lie between 0 and 0.5')
    similarity = np.asarray(similarity)
    if similarity.shape != (len(table), len(table)):
        raise ValueError('similarity must be square and match the table')
    branches = dendrogram(tree, no_plot=True)
    order = np.asarray(branches['leaves'])
    labels = table.id.iloc[order].to_numpy()
    fig = make_subplots(rows=2, cols=2, specs=[[None, {}], [{}, {}]],
                        column_widths=[dendrogram_fraction, 1-dendrogram_fraction],
                        row_heights=[dendrogram_fraction, 1-dendrogram_fraction],
                        horizontal_spacing=0.015, vertical_spacing=0.015)
    # SciPy places leaf centers at 5, 15, ...; map them onto heatmap indices.
    positions, distances = [], []
    for position, distance in zip(branches['icoord'], branches['dcoord']):
        positions.extend([value/10-0.5 for value in position] + [None])
        distances.extend(distance + [None])
    style = dict(mode='lines', line=dict(color='#555555', width=1),
                 showlegend=False, hoverinfo='skip')
    fig.add_trace(go.Scatter(x=positions, y=distances, **style), row=1, col=2)
    fig.add_trace(go.Scatter(x=distances, y=positions, **style), row=2, col=1)
    row_labels = np.broadcast_to(labels[:, None], (len(order), len(order)))
    column_labels = np.broadcast_to(labels[None, :], row_labels.shape)
    fig.add_trace(go.Heatmap(z=similarity[np.ix_(order, order)],
                            x=np.arange(len(order)), y=np.arange(len(order)),
                            customdata=np.stack([row_labels, column_labels], axis=-1),
                            zmin=color_range[0], zmax=color_range[1], colorscale=colorscale,
                            colorbar_title='CKA',
                            hovertemplate='%{customdata[0]}<br>%{customdata[1]}'
                                          '<br>CKA=%{z:.4f}<extra></extra>'), row=2, col=2)
    extent = [-0.5, len(order)-0.5]
    maximum = float(np.max(tree[:, 2]))
    branch_extent = max(maximum*1.05, 1e-6)
    fig.update_xaxes(range=extent, showticklabels=False, row=1, col=2)
    fig.update_yaxes(range=[0, branch_extent], title_text='Linkage distance', row=1, col=2)
    fig.update_xaxes(range=[branch_extent, 0], title_text='Linkage distance', row=2, col=1)
    fig.update_yaxes(range=extent[::-1], showticklabels=False, row=2, col=1)
    fig.update_xaxes(range=extent, tickmode='array', tickvals=np.arange(len(order)),
                     ticktext=labels, showticklabels=show_labels, tickangle=-90, row=2, col=2)
    fig.update_yaxes(range=extent[::-1], tickmode='array', tickvals=np.arange(len(order)),
                     ticktext=labels, showticklabels=show_labels, side='right', row=2, col=2)
    # Match the leaf axes so zooming keeps trees aligned with the heatmap.
    fig.update_xaxes(matches='x3', row=1, col=2)
    fig.update_yaxes(matches='y3', row=2, col=1)
    fig.update_xaxes(showgrid=False, zeroline=False)
    fig.update_yaxes(showgrid=False, zeroline=False)
    fig.update_layout(template='plotly_white', height=height,
                      title='Whole representations: original-space CKA with dendrograms')
    return fig
