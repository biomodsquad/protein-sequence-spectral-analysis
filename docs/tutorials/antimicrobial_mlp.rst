Compare MLPs on antimicrobial peptides
======================================

The existing ABP_GN notebook compares Fourier, circular Z, finite Mellin, and
Haar wavelet magnitudes of corrected BLOSUM100 on Gram-negative antibacterial
peptide classification. It preserves 3,562 supplied training peptides and 1,973
test peptides, with the previous seed-7 stratified training validation split.
The raw CSVs are local inputs and are not redistributed.

.. seealso::

   :download:`Download the editable AMP notebook <../../research/AMP_transform_MLP.ipynb>`
   and read the repository's ``research/README.md`` for execution instructions.

Run the comparison
------------------

Install ``.[ml,wavelets,validation]`` and notebook dependencies. Set ``TRAIN_CSV``
and ``TEST_CSV`` in the controls cell, then run all cells. Defaults look first
in ``research/data/``, then ``~/Downloads/input_csv/``. Required columns are
``Sequence`` and ``Class``; classes are ``Positive`` and ``Negative``.

.. code-block:: bash

   python -m research.execute_amp_notebook \
     --output validation/amp_transform_mlp/my_comparison.executed.ipynb

Choose a fresh output path. The notebook saves fresh timestamped results under
ignored ``validation/amp_transform_mlp/``: metrics, every validation trial,
row-indexed test probabilities, split indices, feature names, source/input hashes,
and environment versions. It plots ROC, precision-recall, and training-loss curves.

What is controlled
------------------

Defaults use 20 channels, a 256-position grid, 5,120 magnitudes per peptide, one
64-unit ReLU MLP candidate, and one model seed. All four models therefore have
327,809 parameters. StandardScaler fits only the fitting partition. Validation
selects candidates; the winner is refit on the full supplied training set. The
MLP's internal early stopping uses training-only data. Threshold stays at 0.5.

This full-coefficient policy differs from the historical dyadic Fourier audit.
Adding candidates or transform settings changes the search budget. Adding model
seeds measures fitting variability, not dataset or family-level uncertainty.
Wavelet ``component='real'`` is a separate experiment because magnitudes lose signs.

Recorded exploratory baseline
-----------------------------

.. list-table:: One seed, existing previously studied test set
   :header-rows: 1

   * - Transform
     - Test balanced accuracy
     - MCC
     - ROC-AUC
     - Average precision
   * - Fourier
     - 76.44%
     - 0.529
     - 0.834
     - 0.749
   * - Z
     - 76.38%
     - 0.525
     - 0.845
     - 0.782
   * - Mellin
     - 74.66%
     - 0.488
     - 0.835
     - 0.787
   * - Wavelet
     - 74.72%
     - 0.493
     - 0.818
     - 0.719

These are recorded runs, not recomputed by building the documentation. The API
unification replay reproduced all test probabilities and metrics. The supplied
test set has already been used in earlier experiments; small score differences
do not establish superiority. A confirmatory study needs repeated seeds and a
family-disjoint or new external test set, with policy fixed before evaluation.
