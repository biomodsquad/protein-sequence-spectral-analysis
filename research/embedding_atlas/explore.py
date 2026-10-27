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
                                              symbol='star' if reference else 'circle')))
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
