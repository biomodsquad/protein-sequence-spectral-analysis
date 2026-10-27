Explore the residue-embedding atlas
===================================

The atlas compares embeddings before sequence-specific transforms. It helps
identify redundant residue-property channels and similar whole representations;
it does not directly measure antimicrobial performance or spectral equivalence.

.. seealso::

   :download:`Download the editable atlas notebook <../../research/embedding_atlas/Embedding_atlas.ipynb>`.
   Work from the checkout so its helpers and reference results are available.

Two views, two meanings
-----------------------

**Individual vectors.** Align channels to the same 20 residues, center and
normalize, and optionally identify sign reversals. The sign-invariant comparison
uses vector outer products, with squared correlation as similarity. Complete
linkage cuts give explicit pairwise similarity guarantees. PCA and UMAP are
visualizations; clustering and neighbors use original-space distances.

**Whole representations.** Compare centered residue Gram matrices with linear
centered kernel alignment, explicitly choosing whether to equalize channel norms.
This handles unequal channel counts. A separate aligned-matrix-score view keeps
signed 20-by-20 substitution/contact scores. Geometry similarity does not imply
equal per-channel Fourier magnitudes: rotating a channel basis can preserve a
Gram matrix while changing channel coefficients.

Inputs and controls
-------------------

Place official AAindex1/2/3 snapshots in ignored ``.aaindex/``. Downloads happen
only after the notebook's explicit ``DOWNLOAD_MISSING=True`` choice. Source
checksum verification is on by default. Raw databases are not distributed.
Install the research requirements documented in
``research/embedding_atlas/README.md``; UMAP and Plotly are not package runtime
requirements or part of the small offline tutorial notebooks.

Edit family/orientation filters, sign conventions, channel weighting, PCA axes,
UMAP neighbors/minimum distance/seed/epochs, clustering linkage/cuts, and display
settings. Refit affected stages after changing inputs. Plot-only edits do not
require a new UMAP fit. Fixed-seed layouts are cached by input/settings.
Exports reject stale fitted stages and preserve a new run's settings/provenance.

Interpretation limits
---------------------

Uniform residue comparisons ignore composition, order, task labels, and
frequency-dependent relevance. Apparent UMAP islands are not quantitative clusters.
Raw-offset equivalence also depends on centering and padding conventions. Treat
atlas-based reduction as a candidate-search heuristic and evaluate the reduced
channels in a training-only pipeline before claiming predictive equivalence.
