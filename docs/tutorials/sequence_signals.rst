From a sequence matrix to a spectrum
====================================

The key idea is a transform **along sequence position**, independently for each
residue-property channel. The amino-acid alphabet is the lookup axis of the
embedding, not the axis being Fourier transformed.

.. figure:: ../_static/figures/sequence-matrix-spectral-analysis.png
   :alt: One-hot rows for I K L F Y R, followed by all-zero padding rows, multiply a 20-by-3 embedding matrix. The resulting three signal columns are Fourier transformed down their length.
   :width: 620px

   Project-provided illustration, reproduced unchanged. Its three channels and
   numerical values illustrate the matrix operation; the code below uses actual
   Kidera factors rather than transcribing the schematic scores.

Read the figure from left to right
----------------------------------

Let a sequence contain :math:`N` residues from an alphabet of size :math:`A=20`.
Its one-hot matrix :math:`O` has shape :math:`N\times20`: each real sequence row
contains one 1 in the column for that residue. An embedding matrix :math:`E`
has shape :math:`20\times C`, with a row of channel scores for each residue.
Multiplication selects the appropriate score row:

.. math::

   X = OE,\qquad X\in\mathbb{R}^{N\times C}.

The figure uses :math:`C=3`. BLOSUM100 supplies 20 channels; Kidera supplies 10.
The embedding can have any supported number of complete residue channels.
Protspect performs the same lookup directly without materializing the one-hot
matrix, which is useful for memory and speed.

Zero-pad the sequence axis
--------------------------

Choose :math:`L\geq N`. Appending :math:`L-N` all-zero rows to the one-hot matrix
produces :math:`\widetilde O\in\mathbb{R}^{L\times20}`. These are padding rows,
not one-hot rows for an extra residue. With raw embedding scores:

.. math::

   \widetilde X = \widetilde O E =
   \begin{bmatrix}X\\0\end{bmatrix}\in\mathbb{R}^{L\times C}.

Padding changes the sampling grid; it does not add residue information. Do not
pad alphabet columns or embedding channels. If you enable centering or channel
normalization, preprocess the **real embedded rows first**, then pad with zeros.
Otherwise padding would contaminate the mean or become nonzero after centering.

Transform each signal column
----------------------------

For the Fourier transform, each channel produces one complex coefficient per
sampled frequency:

.. math::

   F_{k,c} = \sum_{n=0}^{L-1}\widetilde X_{n,c}
       \exp\!\left(-2\pi i\frac{kn}{L}\right).

In array notation this is ``fft(padded_embedding, axis=0)``. ``axis=0`` is the
sequence-length axis, shown by the vertical arrow. ``fftshift(..., axes=0)``
orders coefficients from negative to positive frequency; it does not change
which axis is transformed. Magnitude describes coefficient strength, while phase
retains alignment information that magnitude alone discards.

Check the matrix calculation against Protspect
----------------------------------------------

.. testcode::

   import numpy as np
   from scipy.fft import fft, fftshift
   from protspect import transform_sequence
   from protspect.embeddings import kidera_factors

   sequence = "IKLFYR"
   alphabet = "ACDEFGHIKLMNPQRSTVWY"  # Same column ordering as the figure.
   embedding = dict(list(kidera_factors().items())[:3])
   channels = list(embedding)
   N, L = len(sequence), 16
   one_hot = np.eye(20)[[alphabet.index(aa) for aa in sequence]]
   E = np.array([[embedding[channel][aa] for channel in channels]
                 for aa in alphabet])
   X = one_hot @ E
   padded_one_hot = np.pad(one_hot, ((0, L-N), (0, 0)))
   padded_X = np.pad(X, ((0, L-N), (0, 0)))
   np.testing.assert_allclose(padded_one_hot @ E, padded_X)
   matrix_fft = fftshift(fft(padded_X, axis=0), axes=0)
   result = transform_sequence(sequence, embedding, method="fourier", length=L)
   library_fft = np.column_stack([result.coefficients[c] for c in channels])
   np.testing.assert_allclose(matrix_fft, library_fft)
   assert one_hot.shape == (6, 20) and E.shape == (20, 3)
   assert library_fft.shape == (16, 3)

The same channel signals feed Z, Mellin, and wavelet transforms, although those
backends use different kernels and output coordinates. Continue with
:doc:`representations`, :doc:`transforms`, and :doc:`features`, or open notebook
01 in :doc:`../examples/index` to change the sequence and matrix dimensions.
