"""Compare residue vectors and whole representations without peptide outcome labels.

Run from the repository root with python -m research.embedding_atlas.analyze.
Raw AAindex databases remain user-supplied and are not redistributed.
"""

import argparse
from collections import Counter
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, leaves_list, linkage
from scipy.spatial.distance import squareform
from sklearn.decomposition import PCA
from sklearn.manifold import trustworthiness
from sklearn.neighbors import NearestNeighbors
from threadpoolctl import threadpool_limits

from protspect.aaindex import AAIndex, AAIndex2, AAIndex3
from protspect.embeddings import b100, kidera_factors

ALPHABET = 'ARNDCQEGHILKMFPSTWYV'
COLORS = {'AAindex1': '#4682b4', 'AAindex2': '#e69f00', 'AAindex3': '#009e73',
          'BLOSUM100': 'red', 'Kidera': 'purple'}


# Conventional dominant-property labels; KF5–KF10 mix multiple properties.
# https://search.r-project.org/CRAN/refmans/Peptides/html/kideraFactors.html
KIDERA_DESCRIPTIONS = {
    'KF1': 'Helix/bend preference',
    'KF2': 'Side-chain size',
    'KF3': 'Extended structure preference',
    'KF4': 'Hydrophobicity',
    'KF5': 'Double-bend preference',
    'KF6': 'Partial specific volume',
    'KF7': 'Flat extended preference',
    'KF8': 'Occurrence in alpha region',
    'KF9': 'Carboxyl-group pK (pK-C)',
    'KF10': 'Surrounding hydrophobicity',
}


def unit_profile(values):
    """Remove the equal-residue mean and normalize L2; return None for constants."""
    values = np.asarray(values, dtype=float)
    centered = values-values.mean()
    norm = np.linalg.norm(centered)
    return None if norm <= 1e-12 else centered/norm


def projector_features(profiles):
    """Sign-invariant vector features: vec(v v^T), with unit-Frobenius norm."""
    return np.einsum('ni,nj->nij', profiles, profiles).reshape(len(profiles), -1)


def gram_features(embedding, equal_channels):
    """Centered residue Gram, normalized for linear CKA; never whiten the basis."""
    centered = np.asarray(embedding, dtype=float)-np.mean(embedding, axis=0)
    if equal_channels:
        norms = np.linalg.norm(centered, axis=0)
        centered = centered[:, norms > 1e-12]/norms[norms > 1e-12]
    gram = centered@centered.T
    norm = np.linalg.norm(gram)
    if norm <= 1e-12:
        raise ValueError('representation has no residue contrasts')
    eig = np.linalg.eigvalsh(gram)
    eig = eig[eig > eig.max()*1e-12]
    weights = eig/eig.sum()
    stats = {'rank': len(eig), 'effective_rank': float(np.exp(-np.sum(weights*np.log(weights))))}
    return gram.ravel()/norm, stats


def shape_distance(unit_features):
    """Euclidean distances between unit projectors/Grams; similarity is their dot product."""
    similarity = np.clip(unit_features@unit_features.T, -1, 1)
    distance = np.sqrt(np.maximum(0, 2*(1-similarity)))
    np.fill_diagonal(distance, 0)
    return similarity, distance


def cluster_summary(distance, thresholds, similarity_squared=False):
    """Complete linkage guarantees pairwise similarity at each reported cut."""
    tree = linkage(squareform(distance, checks=False), method='complete', optimal_ordering=False)
    summary, labels = [], {}
    for threshold in thresholds:
        cut = np.sqrt(2*(1-(threshold**2 if similarity_squared else threshold)))
        membership = fcluster(tree, cut, criterion='distance')
        labels[threshold] = membership
        summary.append({'minimum_similarity': threshold, 'clusters': len(set(membership)),
                        'largest_cluster': int(np.bincount(membership).max())})
    return tree, summary, labels


def load_inputs(paths):
    """Keep native axes, exclude only invalid individual vectors or whole matrices."""
    dbs = [cls.from_file(path) for cls, path in zip((AAIndex, AAIndex2, AAIndex3), paths)]
    records, vectors, whole, exclusions = [], [], [], []

    def vector(name, family, description, values, accession, axis='', parent_complete=True):
        if any(v is None for v in values):
            exclusions.append({'kind': 'vector', 'id': name, 'reason': 'missing values'})
            return
        profile = unit_profile(values)
        if profile is None:
            exclusions.append({'kind': 'vector', 'id': name, 'reason': 'constant profile'})
            return
        values = np.asarray(values, dtype=float)
        records.append(dict(id=name, family=family, accession=accession, axis=axis,
                            description=description, mean=values.mean(),
                            population_std=values.std(), parent_complete=parent_complete))
        vectors.append(profile)

    def representation(name, family, description, matrix, orientation, square=False):
        whole.append(dict(id=name, family=family, description=description, embedding=matrix,
                          orientation=orientation, square=square))

    for entry in dbs[0].values():
        vector(entry.accession, 'AAindex1', entry.description,
               [entry.values[a] for a in ALPHABET], entry.accession)
    for family, mapping in [('BLOSUM100', b100), ('Kidera', kidera_factors())]:
        for name, channel in mapping.items():
            description = KIDERA_DESCRIPTIONS[name] if family == 'Kidera' else name
            vector(f'{family}:{name}', family, description,
                   [channel[a] for a in ALPHABET], family)
        e = np.array([[channel[a] for channel in mapping.values()] for a in ALPHABET])
        representation(family, family, family, e, 'builtin', square=family == 'BLOSUM100')
    # An explicitly labelled reference family; normalized-channel CKA prevents raw units dominating.
    complete = [entry for entry in dbs[0].values() if entry.complete]
    representation('AAindex1_complete_set', 'AAindex1', 'All complete scalar indices',
                   np.array([[e.values[a] for e in complete] for a in ALPHABET]), 'all scalars')
    for family, db in zip(('AAindex2', 'AAindex3'), dbs[1:]):
        for entry in db.values():
            if not (set(ALPHABET) <= set(entry.rows) and set(ALPHABET) <= set(entry.columns)):
                exclusions.append({'kind': 'matrix', 'id': entry.accession, 'reason': 'missing alphabet'})
                continue
            m = np.array([[entry.values[r][c] for c in ALPHABET] for r in ALPHABET], dtype=object)
            for axis, array in [('rows', m.T), ('columns', m)]:
                for i, a in enumerate(ALPHABET):
                    vector(f'{entry.accession}:{axis}:{a}', family, entry.description,
                           list(array[:, i]), entry.accession, axis, entry.canonical_complete)
            if entry.canonical_complete:
                m = m.astype(float)
                symmetric = np.array_equal(m, m.T)
                representation(entry.accession+':rows', family, entry.description, m.T,
                               'rows=columns' if symmetric else 'rows', square=True)
                if not symmetric:
                    representation(entry.accession+':columns', family, entry.description, m,
                                   'columns', square=True)
            else:
                exclusions.append({'kind': 'matrix', 'id': entry.accession,
                                   'reason': 'incomplete canonical block'})
    return dbs, pd.DataFrame(records), np.asarray(vectors), whole, pd.DataFrame(exclusions)


def deduplicate(profiles):
    """Collapse sign/affine-equivalent shapes to 10 decimal places; retain all aliases."""
    unique, lookup, assignments = [], {}, []
    for profile in profiles:
        canonical = profile*np.sign(profile[np.argmax(np.abs(profile))])
        key = tuple(np.round(canonical, 10))
        if key not in lookup:
            lookup[key] = len(unique)
            unique.append(profile)
        assignments.append(lookup[key])
    return np.asarray(unique), np.asarray(assignments)


def pca(features):
    """Use exact SVD, with reproducible axes up to standard SVD sign conventions."""
    model = PCA(n_components=2, svd_solver='full')
    coords = model.fit_transform(features)
    return coords, model.explained_variance_ratio_.tolist()


def layout_quality(distance, coords):
    """Full-data neighbor recall plus trustworthiness on a fixed 1,000-point audit subset."""
    k = min(10, len(coords)//3)
    original = NearestNeighbors(n_neighbors=k+1, metric='precomputed').fit(distance)
    hi_raw = original.kneighbors(distance, return_distance=False)
    lo_raw = NearestNeighbors(n_neighbors=k+1).fit(coords).kneighbors(coords, return_distance=False)
    # Explicit self removal remains correct when duplicate whole geometries tie at zero.
    hi = [[j for j in row if j != i][:k] for i, row in enumerate(hi_raw)]
    lo = [[j for j in row if j != i][:k] for i, row in enumerate(lo_raw)]
    recall = np.mean([len(set(a)&set(b))/k for a, b in zip(hi, lo)])
    rows = np.sort(np.random.default_rng(0).choice(len(coords), min(1000, len(coords)), replace=False))
    trust = trustworthiness(distance[np.ix_(rows, rows)], coords[rows], metric='precomputed',
                            n_neighbors=k)
    return {'neighbor_recall_at_10': float(recall), 'audit_trustworthiness_at_10': float(trust),
            'trustworthiness_audit_n': len(rows)}


def scatter(ax, coords, table, title, variance=None):
    """Plot source families and highlight reference representations above dense points."""
    for family, color in COLORS.items():
        rows = np.flatnonzero(table.family.to_numpy() == family)
        if not len(rows):
            continue
        reference = family in ('BLOSUM100', 'Kidera')
        ax.scatter(*coords[rows].T, s=65 if reference else 13, c=color,
                   alpha=1 if reference else .55, marker={'BLOSUM100': '*', 'Kidera': 'X'}.get(family, 'o'),
                   label=f'{family} (n={len(rows)})', edgecolors='black' if reference else 'none',
                   linewidths=0.5 if reference else 0, zorder=4 if reference else 2)
    ax.set_title(title, fontsize=11)
    ax.set_xlabel(f'PC1 ({variance[0]:.1%})' if variance else 'UMAP 1')
    ax.set_ylabel(f'PC2 ({variance[1]:.1%})' if variance else 'UMAP 2')
    ax.grid(alpha=.15)


def save_figure(fig, out, name):
    """Write screen and vector-publication formats."""
    fig.savefig(out/f'{name}.png', dpi=180, bbox_inches='tight')
    fig.savefig(out/f'{name}.pdf', bbox_inches='tight')
    plt.close(fig)


def run(args):
    """Produce complete catalogues, similarities, representatives, layouts, and figures."""
    out = args.output
    out.mkdir(parents=True, exist_ok=False)
    dbs, catalogue, profiles, whole, exclusions = load_inputs(args.inputs)
    print('Loaded', len(catalogue), 'valid vector aliases and', len(whole), 'whole representations', flush=True)
    exclusions.to_csv(out/'excluded.csv', index=False)
    unique, assignments = deduplicate(profiles)
    catalogue['unique_id'] = assignments
    catalogue.to_csv(out/'vector_catalogue.csv', index=False)
    # Prefer reference family colors for shared shapes; aliases always retain every source in CSV.
    priority = {'Kidera': 0, 'BLOSUM100': 1, 'AAindex1': 2, 'AAindex2': 3, 'AAindex3': 4}
    display = (catalogue.assign(priority=catalogue.family.map(priority))
               .sort_values('priority', kind='stable').drop_duplicates('unique_id')
               .set_index('unique_id').loc[np.arange(len(unique))].reset_index())
    display['aliases'] = np.bincount(assignments)
    projectors = projector_features(unique)
    sim, distance = shape_distance(projectors)
    tree, curves, clusters = cluster_summary(distance, [.8, .9, .95, .98, .99], True)
    display['cluster_r095'] = clusters[.95]
    catalogue['cluster_r095'] = clusters[.95][assignments]
    catalogue.to_csv(out/'vector_catalogue.csv', index=False)
    representatives = []
    for cluster in sorted(set(clusters[.95])):
        ids = np.flatnonzero(clusters[.95] == cluster)
        medoid = ids[np.argmin(distance[np.ix_(ids, ids)].mean(axis=1))]
        record = display.iloc[medoid].to_dict()
        record.update(cluster=int(cluster), unique_shapes=len(ids),
                      source_aliases=int(np.bincount(assignments)[ids].sum()))
        representatives.append(record)
    pd.DataFrame(representatives).to_csv(out/'vector_cluster_representatives_r095.csv', index=False)
    a1 = catalogue.family.to_numpy() == 'AAindex1'
    a1_profiles = profiles[a1]
    a1_projectors = projector_features(a1_profiles)
    _, a1_distance = shape_distance(a1_projectors)
    _, a1_curves, a1_clusters = cluster_summary(a1_distance, [.8, .9, .95, .98, .99], True)
    a1_table = catalogue[a1].copy().reset_index(drop=True)
    a1_table['aaindex1_cluster_r095'] = a1_clusters[.95]
    a1_reps = []
    for cluster in sorted(set(a1_clusters[.95])):
        ids = np.flatnonzero(a1_clusters[.95] == cluster)
        medoid = ids[np.argmin(a1_distance[np.ix_(ids, ids)].mean(axis=1))]
        row = a1_table.iloc[medoid].to_dict()
        row['cluster_size'] = len(ids)
        a1_reps.append(row)
    pd.DataFrame(a1_reps).to_csv(out/'aaindex1_representatives_r095.csv', index=False)
    references = catalogue.family.isin(['BLOSUM100', 'Kidera']).to_numpy()
    neighbors = []
    for index in np.flatnonzero(references):
        correlations = a1_profiles@profiles[index]
        for rank, target in enumerate(np.argsort(-abs(correlations))[:3], 1):
            other = a1_table.iloc[target]
            neighbors.append(dict(reference=catalogue.iloc[index].id, family=catalogue.iloc[index].family,
                                  rank=rank, aaindex=other.id, correlation=correlations[target],
                                  absolute_correlation=abs(correlations[target]),
                                  description=other.description))
    pd.DataFrame(neighbors).to_csv(out/'reference_vector_neighbors.csv', index=False)
    vpca, vvar = pca(projectors)
    signed_pca, signed_var = pca(unique)
    display['projector_PC1'], display['projector_PC2'] = vpca.T
    display['signed_PC1'], display['signed_PC2'] = signed_pca.T
    layout_metrics = {'projector_pca': layout_quality(distance, vpca)}
    from umap import UMAP
    layouts = {}
    for nn, seed in [(15, 0), (15, 1), (50, 0)]:
        label = f'umap_n{nn}_seed{seed}'
        print('Fitting', label, 'on', len(unique), 'unique vector shapes', flush=True)
        coords = UMAP(n_components=2, metric='precomputed', n_neighbors=nn, min_dist=.15,
                      n_epochs=300, random_state=seed, n_jobs=1).fit_transform(distance.astype('float32'))
        layouts[label] = coords
        display[label+'_1'], display[label+'_2'] = coords.T
        layout_metrics[label] = layout_quality(distance, coords)
    display.to_csv(out/'unique_vector_layouts.csv', index=False)
    fig, axes = plt.subplots(1, 3, figsize=(18, 6), constrained_layout=True)
    for ax, (name, coordinates) in zip(axes, layouts.items()):
        scatter(ax, coordinates, display, name.replace('_', ' ') +
                f"\ntrustworthiness={layout_metrics[name]['audit_trustworthiness_at_10']:.3f}")
    fig.suptitle('UMAP sensitivity: seed and neighborhood size (same original distances)', fontsize=15)
    save_figure(fig, out, 'umap_sensitivity')
    for name, group in [('vector', curves), ('aaindex1', a1_curves)]:
        for row in group:
            row['scope'] = name
    pd.DataFrame(curves+a1_curves).to_csv(out/'cluster_thresholds.csv', index=False)
    # Whole embeddings: primary equal-channel geometry and secondary native relative weighting.
    grams, raw_grams, rows = [], [], []
    for item in whole:
        gram, stats = gram_features(item['embedding'], True)
        native, native_stats = gram_features(item['embedding'], False)
        grams.append(gram)
        raw_grams.append(native)
        rows.append({k: v for k, v in item.items() if k != 'embedding'} | stats |
                    {'channels': item['embedding'].shape[1],
                     'native_effective_rank': native_stats['effective_rank']})
    grams, raw_grams = np.array(grams), np.array(raw_grams)
    whole_table = pd.DataFrame(rows)
    cka, wd = shape_distance(grams)
    native_cka, _ = shape_distance(raw_grams)
    pd.DataFrame(cka, index=whole_table.id, columns=whole_table.id).to_csv(out/'whole_cka_equal_channels.csv')
    pd.DataFrame(native_cka, index=whole_table.id, columns=whole_table.id).to_csv(out/'whole_cka_native_weights.csv')
    wpca, wvar = pca(grams)
    wumap = UMAP(n_neighbors=15, min_dist=.15, n_epochs=500, metric='precomputed',
                 random_state=0, n_jobs=1).fit_transform(wd.astype('float32'))
    whole_table['PC1'], whole_table['PC2'] = wpca.T
    whole_table['UMAP1'], whole_table['UMAP2'] = wumap.T
    wt, wcurves, wclusters = cluster_summary(wd, [.8, .9, .95, .98, .99])
    whole_table['cluster_cka095'] = wclusters[.95]
    whole_table.to_csv(out/'whole_representations.csv', index=False)
    bindex = whole_table.index[whole_table.id == 'BLOSUM100'][0]
    neighbors = whole_table.copy()
    neighbors['cka_to_b100'] = cka[bindex]
    neighbors['native_cka_to_b100'] = native_cka[bindex]
    neighbors.sort_values('cka_to_b100', ascending=False).to_csv(out/'whole_neighbors_b100.csv', index=False)
    square_rows = np.flatnonzero(whole_table.square)
    aligned = np.array([unit_profile(whole[i]['embedding'].ravel()) for i in square_rows])
    spca, svar = pca(aligned)
    score_correlations = aligned@aligned.T
    pd.DataFrame(score_correlations, index=whole_table.iloc[square_rows].id,
                 columns=whole_table.iloc[square_rows].id).to_csv(out/'aligned_matrix_score_correlations.csv')
    whole_table.iloc[square_rows].assign(score_PC1=spca[:, 0],score_PC2=spca[:, 1]).to_csv(
        out/'aligned_matrix_score_pca.csv', index=False)
    # Plot 1: whole representations and geometry.
    fig, axes = plt.subplots(2, 2, figsize=(15, 12), constrained_layout=True)
    scatter(axes[0, 0], spca, whole_table.iloc[square_rows], 'Aligned matrix scores: signed PCA', svar)
    scatter(axes[0, 1], wpca, whole_table, 'Whole residue geometry: CKA-space PCA', wvar)
    for label in ['BLOSUM100', 'Kidera']:
        i = whole_table.index[whole_table.id == label][0]
        axes[0, 1].annotate(label, wpca[i], xytext=(7, 7), textcoords='offset points', fontsize=10)
    scatter(axes[1, 0], wumap, whole_table, 'Whole residue geometry: UMAP (illustrative)')
    top = neighbors[~neighbors.id.isin(['BLOSUM100','AAindex1_complete_set'])].nlargest(12, 'cka_to_b100')
    axes[1, 1].barh(top.id.iloc[::-1], top.cka_to_b100.iloc[::-1], color='#4682b4')
    axes[1, 1].set_xlim(0, 1)
    axes[1, 1].set_xlabel('Equal-channel linear CKA with BLOSUM100')
    axes[1, 1].set_title('Closest whole representations to BLOSUM100')
    axes[0, 0].legend(fontsize=8, loc='best')
    axes[0, 1].legend(fontsize=8, loc='best')
    fig.suptitle('Whole representations across the same 20 canonical residues', fontsize=16)
    save_figure(fig, out, 'whole_overview')
    # Plot 2: individual vectors and sign-invariant projectors.
    fig, axes = plt.subplots(2, 2, figsize=(15, 12), constrained_layout=True)
    scatter(axes[0, 0], signed_pca, display, 'Signed-profile PCA (one sign representative per shape)', signed_var)
    scatter(axes[0, 1], vpca, display, 'Sign-invariant projector PCA', vvar)
    scatter(axes[1, 0], layouts['umap_n15_seed0'], display, 'Sign-invariant UMAP (n_neighbors=15, seed=0)')
    for curve, label in [(curves, 'All unique channel shapes'), (a1_curves, 'Complete AAindex1 channels')]:
        axes[1, 1].plot([r['minimum_similarity'] for r in curve], [r['clusters'] for r in curve],
                        'o-', label=label)
    axes[1, 1].set(xlabel='Minimum within-cluster |Pearson r|', ylabel='Complete-linkage clusters',
                    title='Conservative candidate-reduction thresholds')
    axes[1, 1].legend(fontsize=9)
    axes[0, 0].legend(fontsize=8)
    fig.suptitle('Individual channels: centered shape comparison, sign-equivalent for amplitudes', fontsize=15)
    save_figure(fig, out, 'vectors_overview')
    fig, axes = plt.subplots(1, 2, figsize=(16, 7), constrained_layout=True)
    order = leaves_list(wt)
    for ax, matrix, title in [(axes[0], cka, 'Equal-channel weighting'),
                               (axes[1], native_cka, 'Native relative channel weighting')]:
        im = ax.imshow(matrix[np.ix_(order, order)], vmin=0, vmax=1, cmap='viridis', rasterized=True)
        ax.set(title=title, xlabel='Representations, same cluster order', ylabel='Representations')
        ax.set_xticks([])
        ax.set_yticks([])
        pos = int(np.flatnonzero(order == bindex)[0])
        ax.axhline(pos, color='#ff66aa', lw=.8)
        ax.axvline(pos, color='#ff66aa', lw=.8)
    fig.colorbar(im, ax=axes, label='Linear CKA', shrink=.8)
    fig.suptitle('Whole-representation similarity; pink cross marks BLOSUM100', fontsize=15)
    save_figure(fig, out, 'whole_similarity_heatmaps')
    ref = pd.read_csv(out/'reference_vector_neighbors.csv')
    ref = ref[ref['rank'] == 1]
    fig, ax = plt.subplots(figsize=(12, 10), constrained_layout=True)
    labels = [f'{r.reference}  →  {r.aaindex}' for r in ref.itertuples()]
    ax.barh(labels[::-1], ref.absolute_correlation.iloc[::-1],
            color=[COLORS[f] for f in ref.family.iloc[::-1]])
    ax.set(xlim=(0, 1), xlabel='Absolute Pearson correlation across 20 residues',
           title='Closest complete AAindex1 property to each BLOSUM100 / Kidera channel')
    ax.axvline(.95, ls='--', color='gray', lw=1)
    save_figure(fig, out, 'reference_channel_neighbors')
    # Local caches stay outside commits; no source tables are bundled in distributions.
    np.savez_compressed(out/'pairwise_local.npz', vector_distance=distance,
                        whole_cka=cka, whole_native_cka=native_cka)
    manifest = {'inputs': [{'section':i+1,'path':str(Path(p).resolve()),'sha256':d.sha256,'records':len(d)}
                           for i,(p,d) in enumerate(zip(args.inputs,dbs))],
                'analysis_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                'builtin_sha256': hashlib.sha256(json.dumps({'b100':{k:dict(v) for k,v in b100.items()},
                                    'kidera':kidera_factors()},sort_keys=True).encode()).hexdigest(),
                'vector_aliases':len(catalogue),'unique_shapes':len(unique),
                'vector_families':dict(Counter(catalogue.family)),
                'vector_clusters':curves,'aaindex1_clusters':a1_curves,
                'whole_representations':len(whole),'whole_clusters':wcurves,
                'vector_pca_variance':vvar,'signed_pca_variance':signed_var,
                'whole_pca_variance':wvar,'score_pca_variance':svar,
                'layout_metrics':layout_metrics,
                'whole_layout_metrics':layout_quality(wd,wumap),
                'environment':{'python':platform.python_version(), **{p:importlib.metadata.version(p)
                                for p in ['numpy','scipy','pandas','matplotlib','scikit-learn','umap-learn']}},
                'parameters':{'residue_order':ALPHABET,'dedup_decimals':10,'vector_threshold':.95,
                              'umap_epochs':300,'umap_min_dist':.15,'threads':1}}
    (out/'report.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps({k:manifest[k] for k in ['vector_aliases','unique_shapes','whole_representations',
                                             'aaindex1_clusters','layout_metrics']},indent=2),flush=True)
    return manifest


def main():
    """Require local sources and a fresh destination; all sources remain offline."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs', nargs=3, type=Path, required=True, metavar=('AAINDEX1','AAINDEX2','AAINDEX3'))
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    with threadpool_limits(limits=1):
        run(args)


if __name__ == '__main__':
    main()
