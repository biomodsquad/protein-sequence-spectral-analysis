Export and select coefficient features
======================================

The shared exporter makes information choices explicit. ``amplitude`` is an
absolute magnitude; ``phase`` is an angle in radians; ``real`` and ``imaginary``
retain signed components. Phase at zero amplitude has no useful interpretation.
For real wavelets, magnitude loses signs and phase is mostly a sign encoding.

.. testcode::

   import numpy as np
   from protspect import transform_sequence

   result = transform_sequence("IKLFYR", method="fourier", length=16)
   values, names = result.features(("amplitude", "phase"),
       selector=lambda c: c["frequency"] >= 0)
   assert values.shape == (20 * 8 * 2,)
   assert names.shape == values.shape
   # Original coefficient indices remain in names after selection.
   assert names[0].endswith(":8")

The selector receives a dictionary, not an unlabelled frequency vector. The
returned boolean mask must match coefficient width. Wavelet selection instead
uses ``lambda c: c['band'] == 'D1'``. Mellin selection uses ``log_frequency``.
A frequency selector is meaningful only on an actual frequency axis:

.. testcode::

   from protspect.utility import frequency_select_physical

   values, names = result.features(selector=lambda c: frequency_select_physical(c["frequency"]))
   assert values.shape == (20 * 4,)  # DC, 1/16, 1/8, 1/4; no repeated +1/2 bin.

Within-sequence channel products
--------------------------------

.. testcode::

   cross = result.cross()
   np.testing.assert_allclose(cross.coefficients["A_R"],
       result.coefficients["A"] * result.coefficients["R"].conj())
   features, names = result.features(kinds=("coefficients", "cross_coefficients"))
   assert features.shape == ((20 + 210) * 16,)

Within-result products use the upper channel triangle, including auto-products.
Their magnitudes obey :math:`|X_i\overline{X_j}|=|X_i||X_j|`, so amplitude-only
products add magnitude interactions but discard relative phase. Fourier products
are not averaged cross-spectral-density estimates. Wavelet products match bands
and positions; Mellin products match log-frequency coordinates.

Collections and metadata
------------------------

``SequenceSet.export_features`` accepts the same component, kind, and selector
controls. ``attributes`` adds metadata columns with collision checks. Every item
shares a grid and schema. For sklearn use :doc:`machine_learning` so the selector
is evaluated once at fit rather than independently for each sample. Feature names
encode channel, method, component, and original index; use the stored coordinates
for the actual frequency or band meaning.

Primary exports contain all unique coefficients unless selected explicitly.
Historical ``export_datatable`` keeps its closed Fourier grid and ``Amp_``/``Phase_``
names. Read :doc:`../UNIFIED_API` before replacing columns for an existing model.
