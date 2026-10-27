# Unsupervised residue-embedding atlas

Compare individual residue-property vectors and whole representations across
AAindex1/2/3, the ten Kidera factors, and Protspect's corrected BLOSUM100.
No peptide labels, predictive-model fitting, or performance claims enter this
analysis. The committed [results](../../validation/embedding_atlas/RESULTS.md)
include quantitative similarities, cluster representatives, PCA/UMAP layouts,
and PNG/PDF figures. Raw AAindex files are not redistributed.

## Explore in a notebook

Open [Embedding_atlas.ipynb](Embedding_atlas.ipynb) in Jupyter with this checkout
as the working directory. Install the research-only notebook dependencies:

```sh
python -m pip install -r research/embedding_atlas/notebook_requirements.txt
```

Run cells in order once. Edit the configuration cells to choose source families,
row/column orientations, signed versus sign-invariant vector comparisons,
deduplication, whole-representation channel weighting, PCA components/axes,
UMAP neighbors/minimum distance/seed/epochs, and clustering method/cuts.
Colors, marker sizes, opacity, heatmap scales, and reference lookups are editable.
Vector centering and unit normalization are fixed comparison conventions.

Rerun the affected stage and its downstream cells after changing inputs. Plot-only
changes do not require refitting. Fixed-seed UMAP layouts are cached in memory by
input distances and effective parameters; restarting the kernel clears that cache.
Cluster cuts and nearest neighbors always use the original geometry.

The notebook reads `.aaindex/aaindex1`, `aaindex2`, and `aaindex3`. Set
`DOWNLOAD_MISSING=True` to explicitly fetch absent snapshots. Reference checksum
verification defaults to enabled; disable it deliberately when exploring a newer
snapshot. Raw downloads, kernel caches, and optional exported runs stay under the
ignored `.aaindex/` directory. Set `EXPORT_RUN=True` in the last cell to save a fresh
run with coordinates, clusters, medoids, similarities, controls, and source hashes.
The export cell rejects stale PCA, UMAP, or cluster results after input changes.

For a fully executed local copy, including plots:

```sh
python -m research.embedding_atlas.execute_notebook \
  --output validation/embedding_atlas/Embedding_atlas.executed.ipynb
```

Choose a fresh output filename for later executions. The large executed notebook
is a local artifact; the clean editable notebook and execution helper are tracked.

## Reproduce

Use an isolated Python 3.12 environment and install this branch plus the analysis
requirements. The listed top-level/numerical dependencies are pinned; this is
not a complete lock of every transitive package. Full input checksums and runtime
versions are recorded in `report.json`.

```sh
python -m venv .venv-atlas
.venv-atlas/bin/python -m pip install -r research/embedding_atlas/requirements.txt
.venv-atlas/bin/python -m pip install --no-deps -e .
.venv-atlas/bin/python -m pytest -q research/embedding_atlas/test_analysis.py
```

Download official AAindex snapshots explicitly using the package helpers, then
run with a fresh output directory (never overwrite an earlier result):

```sh
MPLCONFIGDIR=/private/tmp/protspect-atlas-mpl \
NUMBA_CACHE_DIR=/private/tmp/protspect-atlas-numba \
NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
.venv-atlas/bin/python -m research.embedding_atlas.analyze \
  --inputs .aaindex/aaindex1 .aaindex/aaindex2 .aaindex/aaindex3 \
  --output validation/embedding_atlas_replay

.venv-atlas/bin/python -m research.embedding_atlas.interactive \
  validation/embedding_atlas_replay
```

Use an appropriate writable cache path on other operating systems. The base
package does not depend on UMAP, Plotly, or this research environment. The optional
HTML viewer bundles Plotly for offline hover, zoom, identifier lookup, alias
inspection, and whole-representation neighbors. Its large generated file and
the pairwise cache are local artifacts; regenerate the viewer from committed CSVs.

## Individual vectors

1. Align every channel to `ARNDCQEGHILKMFPSTWYV` with equal residue weights.
   Include AAindex1, Kidera, BLOSUM100, and both row and column profiles of matrix
   entries, including any individually complete profiles from otherwise incomplete
   matrices. Exclude missing or constant individual profiles, recording reasons.
2. Center each vector across the 20 residue types, then normalize its L2 norm.
   Preserve original mean and standard deviation in the catalogue.
3. Deduplicate normalized shapes up to sign, rounded to ten decimal places;
   keep every accession/axis alias. Unique points prevent exact duplicates from
   dominating layouts. For a shared point, displayed source color prioritizes
   Kidera, BLOSUM100, AAindex1, AAindex2, then AAindex3. It is not an exclusive
   source classification; `vector_catalogue.csv` retains all aliases.
4. Use unit outer products `P=v vᵀ` as sign-invariant comparison features.
   Their dot product is `r²`, and their Euclidean distance is `sqrt(2(1-r²))`.
   This produces coherent distances for PCA, UMAP, and hierarchical clustering
   without arbitrary sign flips at projection boundaries. An additional signed
   profile PCA retains one original sign representative per deduplicated shape;
   it is illustrative and depends on that representative choice.
5. Complete-linkage clustering at a correlation threshold guarantees every pair
   in a cluster has at least that absolute correlation (within numerical precision).
   Export distance medoids, cluster sizes, and alias mappings. AAindex1-only
   clustering is performed separately, so its reduction does not depend on how
   many matrix channels were added.

Sign reversal alone preserves amplitudes. Centering and unit scaling, however,
are **comparison choices**, not declarations of identical raw spectral features.
A constant offset can affect nonzero-frequency bins when a finite sequence is
zero-padded. Affine-related channels become equivalent in amplitude after
per-sequence centering plus per-channel energy normalization, apart from
zero-energy cases. Phase and signed cross-spectrum interpretations differ.

## Whole representations

Canonical matrix blocks supply rows as channels; distinct column orientations
are included for asymmetric matrices. A symmetric matrix is represented only
once. Omit incomplete canonical blocks. Include BLOSUM100, Kidera, and the
complete AAindex1 collection as reference representations.

Two complementary comparisons are provided:

- **Aligned matrix-score profiles:** center and unit-normalize the flattened
  canonical matrix; preserve sign and row/column labels. PCA and signed score
  correlations apply to square matrices with matching 20-channel reference axes.
  Kidera's ten channels and the complete AAindex1 family are excluded from this
  score-alignment view, because their channel axes have different meanings.
- **Residue geometry:** construct `E` with 20 residue rows and channel columns,
  center columns across residues, and compare unit-Frobenius residue Gram
  matrices `K=E Eᵀ`. Their cosine similarity is linear centered kernel alignment
  (CKA). This accommodates different channel counts. The primary analysis
  equalizes each channel's norm before forming K; a sensitivity analysis keeps
  native relative channel weights. CKA's orthogonal-channel-rotation invariance
  applies to E after whichever preprocessing was chosen. No whitening is used.

CKA similarity does not imply identical per-channel Fourier amplitudes: changing
basis can preserve the Gram matrix while changing channel magnitudes. Use the
individual-vector comparison alongside whole geometry. Full-span principal-angle
comparisons were avoided because many full-rank centered representations span
the same 19-dimensional space and would appear indistinguishable. Report matrix
rank and entropy-based effective rank instead.

## Layout diagnostics and limitations

Clusters and nearest neighbors use original distances, never distances in a
2D layout. Exact-SVD PCA reports explained variance. UMAP uses precomputed
projector distances with 300 epochs, `min_dist=.15`, seeds 0/1 at 15 neighbors,
and seed 0 at 50 neighbors. Whole-geometry UMAP uses seed 0, 15 neighbors, and
500 epochs. Compute full-data 10-neighbor recall and trustworthiness on a fixed
1,000-point subset (all points for whole representations). Repeated seeds and
neighborhood sizes assess visualization sensitivity, not scientific uncertainty.
UMAP island size, spacing, and apparent boundaries are not quantitative clusters.

These comparisons use uniform residue weights and selected centering/scale
invariances. They do not assess dataset-specific composition, sequence order,
frequency information, supervised relevance, or generalization. Two-dimensional
plots retain only a fraction of the original geometry. Candidate reduction
remains a heuristic to benchmark, not evidence of interchangeable predictors.

Methods reference: Kornblith et al. (2019),
[Similarity of Neural Network Representations Revisited](https://proceedings.mlr.press/v97/kornblith19a.html).
See [AAindex attribution](../../docs/AAINDEX.md#database-attribution) and
[Kidera provenance](../../docs/REPRESENTATIONS.md).
