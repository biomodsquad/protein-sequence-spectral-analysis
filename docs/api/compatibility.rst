Fourier compatibility adapters
==============================

Existing calls preserve Fourier values, closure bins, selector behavior, and
``Amp_``/``Phase_`` names. New code should follow :doc:`../UNIFIED_API` to select the
primary full-bin schema deliberately. Recompute and refit before changing a
trained model's feature definition.

.. autoclass:: protspect.features.SequenceFFT
   :members: compute_spectra, compute_cross_spectra, compute_spectral_properties
   :inherited-members:
   :exclude-members: compute_transform, features, embed_sequence

.. autoclass:: protspect.sklearn.SpectralTransformer

.. automethod:: protspect.features.SequenceSet.export_datatable

.. automethod:: protspect.features.ComplexSet.export_datatable
