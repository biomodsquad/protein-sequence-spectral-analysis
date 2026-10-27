SequenceSet and ComplexSet
==========================

Both collection classes support every transform. ``from_sequence_list`` stores
results and returns None; ``export_features`` is the shared export method.
``spectrum_length`` and ``export_datatable`` remain compatibility spellings.

.. autoclass:: protspect.features.SequenceSet
   :members: from_sequence_list, export_features

.. autoclass:: protspect.features.ComplexSet
   :members: from_sequence_list, export_features
