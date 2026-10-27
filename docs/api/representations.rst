Representations and selectors
=============================

Built-in channels
-----------------

``protspect.embeddings.b100`` is the read-only corrected 20-channel BLOSUM100
mapping. ``protspect.embeddings.kidera`` is the read-only ten-factor mapping.
Use the factory below for independent dictionaries, including sklearn cloning.

.. autofunction:: protspect.embeddings.kidera_factors

AAindex catalogues
------------------

No catalogue is downloaded at import time. Keep local inputs and hashes.

.. autoclass:: protspect.aaindex.AAIndex
   :members: from_file, search, embedding, metadata

.. autoclass:: protspect.aaindex.AAIndexEntry
   :members: values, missing_residues, complete

.. autoclass:: protspect.aaindex.AAIndex2
   :members:
   :inherited-members:
   :exclude-members: get, items, keys, values

.. autoclass:: protspect.aaindex.AAIndex3
   :members:
   :inherited-members:
   :exclude-members: get, items, keys, values

.. autoclass:: protspect.aaindex.AAIndexMatrixEntry
   :members: values, shape, missing_pairs, complete, canonical_complete

.. autofunction:: protspect.aaindex.download_aaindex1

.. autofunction:: protspect.aaindex.download_aaindex2

.. autofunction:: protspect.aaindex.download_aaindex3

Frequency helpers
-----------------

Wrap these in a coordinate-dictionary selector for the primary API, for example
``selector=lambda c: frequency_select_positive(c['frequency'])``. They are not
wavelet band or Mellin log-frequency selectors.

.. autofunction:: protspect.utility.frequency_select_physical

.. autofunction:: protspect.utility.frequency_select_positive
