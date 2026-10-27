Getting started
===============

Install the development API
---------------------------

Use Python 3.10 or newer. The published alpha retains the older Fourier API;
these tutorials use the unified development code. Clone the repository and
install into a virtual environment:

.. code-block:: bash

   git clone https://github.com/biomodsquad/protein-sequence-spectral-analysis.git
   cd protein-sequence-spectral-analysis
   python -m venv .venv
   # macOS/Linux:
   source .venv/bin/activate
   # Windows PowerShell instead: .venv\Scripts\Activate.ps1
   python -m pip install -e '.[ml,wavelets,validation]'

Before the development changes are merged to main, check out the unified API /
documentation branch or its corresponding pull request before installing.
The core package needs NumPy, SciPy, and pandas. The ``ml`` extra adds sklearn;
``wavelets`` adds PyWavelets; ``validation`` adds plotting and research utilities.
There are no implicit AAindex downloads.

For notebooks also install:

.. code-block:: bash

   python -m pip install jupyterlab nbformat nbclient ipykernel
   python -m jupyter lab

For the published Fourier-only alpha, explicitly use
``python -m pip install Protspect==0.1.0a1`` and its existing compatibility calls.
Do not expect the development examples to run against that release.

Make your first table
---------------------

.. testcode::

   from protspect import SequenceSet

   dataset = SequenceSet(method="fourier", compute_cross_spectra=False)
   dataset.from_sequence_list(["ACDE", "ACDF"], length=8,
       attributes={"id": ["peptide-1", "peptide-2"]})
   table = dataset.export_features("amplitude", attributes=["id"])
   assert table.shape == (2, 161)  # ID + 20 channels × 8 bins
   assert list(table.id) == ["peptide-1", "peptide-2"]

Inputs are case-sensitive. Unknown residues, empty sequences, mismatched pair
counts, and grids shorter than a chain raise errors. Clean and validate your
sequence alphabet deliberately rather than silently replacing unknown residues.

Switch transforms
-----------------

.. testcode::

   from protspect import transform_sequence

   wavelet = transform_sequence("ACDE", method="wavelet", length=16,
       wavelet="haar", level=2, mode="periodization")
   values, names = wavelet.features("real")
   assert values.shape == (320,)
   assert set(wavelet.coordinates["band"]) == {"A2", "D2", "D1"}

For sets and sklearn, backend options go in ``transform_options``. Components,
selection, and product export use the same conventions for all methods.

Choose a learning path
----------------------

* **First feature workflow:** :doc:`framework`, :doc:`tutorials/sequence_signals`,
  and notebook 01 in :doc:`examples/index`.
* **Compare transforms or paired chains:** :doc:`tutorials/transforms`,
  :doc:`tutorials/paired_chains`, and notebook 02.
* **Train models:** :doc:`tutorials/machine_learning`, notebook 03, then
  :doc:`tutorials/antimicrobial_mlp` with the existing ABP_GN files.
* **Expand representations:** :doc:`tutorials/representations` and
  :doc:`tutorials/embedding_atlas`.
* **Update older code:** :doc:`UNIFIED_API` and :doc:`COMPATIBILITY`.
