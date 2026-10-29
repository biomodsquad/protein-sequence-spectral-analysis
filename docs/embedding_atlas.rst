Embedding Atlas
===============

Explore how amino-acid representations relate **before applying a sequence
transform**. This atlas compares individual property channels and the complete
residue geometry of embedding matrices. It includes BLOSUM100, the ten Kidera
factors, and valid AAindex1/2/3 properties and matrix profiles.

Hover to inspect identifiers and descriptions, drag to zoom, double-click to
reset, and click legend entries to hide or show families. BLOSUM100 appears as
**red stars** and Kidera as **purple X markers**, both with thin black outlines.
The figures below are exports of the current :doc:`tutorials/embedding_atlas`
notebook at its default settings; changing hyperparameters requires rerunning
that notebook. No code cells are included here.

How the analysis was performed
------------------------------

* Each individual channel is centered across the 20 canonical amino acids and
  scaled to unit L2 norm. Default sign-invariant comparison uses squared Pearson
  correlation, so opposite-sign profiles have the same shape. Duplicate shapes
  share a point; hover text retains their aliases. This is normalization over
  amino acids, with no peptide-data fitting or outcome labels.
* For whole representations, each channel is centered and normalized separately
  to give it equal weight. The residue Gram matrix is then Frobenius-normalized;
  its inner product with another normalized Gram gives linear **CKA**. This
  compares residue geometry even when channel counts differ; it does not whiten
  the embedding. Constant channels do not contribute.
* PCA is fitted to vector projectors or normalized Gram features. UMAP uses
  their original-space distances, ``sqrt(2 * (1 - similarity))``, rather than PCA
  coordinates. UMAP uses 15 neighbors, minimum distance 0.15, seed 0, and 300
  epochs for channels or 500 for whole representations. PCA fits ten components
  and displays the first two; axis labels report explained variance.
* Hierarchical clustering uses complete linkage on the original-space distances.
  Similarity, clustering, and heatmap values are independent of the 2D layouts.
  Distances on PCA/UMAP plots can distort the full geometry and do not measure
  predictive performance.

Individual property channels
----------------------------

PCA
~~~

Each point is a unique centered property shape across amino acids. PCA uses
sign-invariant projector features, so properties with strong positive or negative
correlation are close in the comparison space. Family labels identify a preferred
source for shared shapes, rather than exclusive ownership of that property.

.. raw:: html

   <iframe src="_static/embedding-atlas/vector-pca.html" title="Individual-channel PCA" loading="lazy" style="width:100%;height:680px;border:0;"></iframe>

UMAP
~~~~

UMAP emphasizes local neighborhoods in the same sign-invariant original space.
Hover on Kidera points to see the property descriptions, such as hydrophobicity,
alongside stable factor identifiers. Islands and visual gaps depend on the UMAP
settings; use original-space similarities for quantitative comparisons.

.. raw:: html

   <iframe src="_static/embedding-atlas/vector-umap.html" title="Individual-channel UMAP" loading="lazy" style="width:100%;height:680px;border:0;"></iframe>

Whole representations and matrix geometry
-----------------------------------------

CKA-space PCA
~~~~~~~~~~~~~

Each point now represents a complete embedding configuration. AAindex2/3 matrices
supply canonical row or column profiles; asymmetric orientations can therefore
appear separately. Kidera has ten channels, BLOSUM100 has twenty, and the complete
AAindex1 property collection is another configuration. CKA compares the residue
relationships encoded by these different channel sets.

.. raw:: html

   <iframe src="_static/embedding-atlas/whole-pca.html" title="Whole-representation CKA-space PCA" loading="lazy" style="width:100%;height:680px;border:0;"></iframe>

CKA-space UMAP
~~~~~~~~~~~~~~

The same configurations are placed using the original CKA-derived distances.
This view emphasizes local geometric similarity without requiring equal numbers
of channels or a shared channel basis.

.. raw:: html

   <iframe src="_static/embedding-atlas/whole-umap.html" title="Whole-representation CKA-space UMAP" loading="lazy" style="width:100%;height:680px;border:0;"></iframe>

Original-space CKA heatmap and dendrograms
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Colors show the actual pairwise CKA values, not distances inferred from a
projection. Both axes follow the same complete-linkage leaf order; the top and
left dendrograms show the corresponding original-space linkage distances.
Hover to inspect both representation identifiers and their similarity. Zooming
the heatmap keeps the corresponding dendrogram leaf axes aligned.

.. raw:: html

   <iframe src="_static/embedding-atlas/cka-heatmap.html" title="Original-space CKA heatmap with row and column dendrograms" loading="lazy" style="width:100%;height:960px;border:0;"></iframe>

Inspect the embedding matrix
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This default notebook view shows BLOSUM100 with sequence residues as rows and
reference channels as columns. The colors are the **raw signed scores**, without
channel normalization. This inspection panel is distinct from the normalized
geometry comparisons above; the notebook allows selecting other representations.

.. raw:: html

   <iframe src="_static/embedding-atlas/channel-matrix.html" title="BLOSUM100 raw embedding channel matrix" loading="lazy" style="width:100%;height:610px;border:0;"></iframe>

Signed aligned matrix-score PCA
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Only square matrices with matching canonical residue axes enter this comparison.
Each matrix is flattened after global mean removal and unit-norm scaling. This
preserves channel labels, signs, and relative score patterns, unlike the CKA
comparison of residue geometry. Kidera and the complete scalar-index collection
are excluded because their channel axes do not match the square matrices.

.. raw:: html

   <iframe src="_static/embedding-atlas/matrix-pca.html" title="Signed aligned matrix-score PCA" loading="lazy" style="width:100%;height:680px;border:0;"></iframe>

Sources and reproducibility
---------------------------

The AAindex files are explicitly downloaded from the official database and
SHA256-checked before analysis. Raw database files are not bundled in this site.
These are derived exploratory plots, not a supervised feature-selection result.
Kidera labels describe dominant properties; factors 5–10 combine multiple
properties. See :doc:`REPRESENTATIONS`, :doc:`AAINDEX`, and
:doc:`AAINDEX_MATRICES` for sources and interpretation.

The :doc:`tutorials/embedding_atlas` page links the editable notebook and explains
its controls. The :download:`plot export manifest <_static/embedding-atlas/manifest.json>`
records source, data-snapshot, and figure hashes. To explore sequence spectral
features after choosing channels, continue with :doc:`tutorials/transforms`.
