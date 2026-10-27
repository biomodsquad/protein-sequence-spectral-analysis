Examples and notebooks
======================

The first three notebooks are self-contained offline tutorials. Download a
writable copy and use a Python kernel with the development package plus
``.[ml,wavelets,validation]`` installed. They require no external peptide dataset,
AAindex download, UMAP, or MISTIC. Their code is executed in documentation CI.

After installation, locate the bundled notebooks without assuming a
site-packages path, then copy one to a writable workspace:

.. testcode::

   from importlib.resources import files
   import shutil

   examples = files("protspect.examples")
   notebook = examples.joinpath("01_matrix_to_features.ipynb")
   assert notebook.is_file()
   # shutil.copyfile(notebook, "my_first_tutorial.ipynb")

1. Matrix view and first features
---------------------------------

:download:`01 — Matrix view and first features <../../protspect/examples/01_matrix_to_features.ipynb>`
   Read the supplied one-hot/embedding diagram, build the matrices explicitly,
   check the axis-0 FFT against Protspect, compare BLOSUM100 and Kidera, and
   select coefficients from a sequence table. Edit the sequence and padded length.

2. Transforms and paired chains
-------------------------------

:download:`02 — Transforms and paired chains <../../protspect/examples/02_transforms_and_pairs.ipynb>`
   Compare four transforms of one signal, adjust Z/Mellin/wavelet settings,
   inspect actual axes, and use ComplexSet for aligned cross-chain products.
   Includes an executable conjugation check when chains are reversed.

3. Training-only machine learning
---------------------------------

:download:`03 — Fitted pipelines and synthetic comparison <../../protspect/examples/03_fitted_pipeline.ipynb>`
   Generate labelled synthetic sequences, split before fitting, tune four
   transforms in a Pipeline, and evaluate the selected model on held-out data.
   Inspect frozen channel moments and JSON-safe metadata. Synthetic accuracy is
   not evidence of antimicrobial or other biological performance.

Research workflows
------------------

These notebooks use checkout helpers or user-managed inputs. They are not run
as part of the small offline documentation suite.

:download:`AMP four-transform MLP benchmark <../../research/AMP_transform_MLP.ipynb>`
   Existing ABP_GN split, matched baseline MLP sizes, validation selection,
   training-only scaling, ROC/PR/loss curves, and provenance exports.
   Start with :doc:`../tutorials/antimicrobial_mlp`.

:download:`Embedding atlas <../../research/embedding_atlas/Embedding_atlas.ipynb>`
   Compare AAindex, BLOSUM100, and Kidera channels and whole representations.
   Requires explicit local AAindex snapshots and the research plotting stack.
   Start with :doc:`../tutorials/embedding_atlas`.

:download:`Existing transform sandbox <../../examples/alternative_transforms.ipynb>`
   A compact plotting notebook for transform settings and fixed-width features.

Small Python examples
---------------------

Run these as modules from a checkout, or with the editable package installed.
The AAindex examples require explicit local inputs.

.. literalinclude:: ../../examples/basic_features.py
   :language: python
   :caption: examples/basic_features.py

:download:`Basic collection example <../../examples/basic_features.py>` ·
:download:`Normalized sklearn pipeline <../../examples/normalized_pipeline.py>` ·
:download:`AAindex1 adapter <../../examples/aaindex_features.py>` ·
:download:`AAindex2/3 adapters <../../examples/aaindex_matrix_features.py>`

The old CPP/MISTIC research notebook remains under ``research/`` with its own
locked environment and data provenance. It is a separate scientific example,
not an installation prerequisite for Protspect.
