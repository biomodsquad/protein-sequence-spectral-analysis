Paired chains with ComplexSet
=============================

``ComplexSet`` remains the collection interface for pairs and supports every
transform. It uses one grid across both chain collections, preprocesses each
chain independently, and stores chain A, chain B, and ordered paired results.

.. testcode::

   import numpy as np
   from protspect import ComplexSet, transform_pair

   pairs = ComplexSet(method="z", transform_options={"radius": 1.05})
   pairs.from_sequence_list(["IKLFYR", "ACDE"], ["ACDE", "IKLFYR"], length=16,
       attributes={"pair_id": ["p1", "p2"]})
   paired = pairs.export_features(("amplitude", "phase"), attributes=["pair_id"])
   chain_a = pairs.export_features("real", chain_ID="A")
   chain_b = pairs.export_features("real", chain_ID="B")
   assert paired.shape == (2, 400 * 16 * 2 + 1)
   assert chain_a.shape == chain_b.shape == (2, 20 * 16)

The paired coefficient is :math:`A_i\overline{B_j}` for every ordered channel
combination, not only an upper triangle. Swapping chains conjugates the
corresponding swapped-channel product:

.. testcode::

   forward = transform_pair("IKLFYR", "ACDE", method="fourier", length=16)
   backward = transform_pair("ACDE", "IKLFYR", method="fourier", length=16)
   np.testing.assert_allclose(forward.coefficients["A-A_B-R"],
                             backward.coefficients["A-R_B-A"].conj())

These products describe aligned coefficient coordinates. They do not align two
structures or infer residue contacts. Sequence origin, length, padding, and
relative alignment still matter. Wavelet products are band/position descriptors;
they are not Fourier cross-spectra. Amplitude-only paired products lose phase
and chain-order information in the corresponding swapped channel combination.

Use paired samples in a Pipeline
--------------------------------

.. testcode::

   from protspect.sklearn import SequenceTransformer

   features = SequenceTransformer(method="wavelet", length=16, paired=True,
       component="real", transform_options={"wavelet": "haar", "level": 2})
   X = features.fit_transform([("IKLFYR", "ACDE"), ("ACDE", "IKLFYR")])
   assert X.shape == (2, 400 * 16)

With ``paired=True`` each sample must contain exactly two sequence strings.
Fitted channel moments, if requested, weight every real residue from both
training chains equally. Keep related pairs and repeated chains grouped during
splitting when your scientific question requires independence at that level.
