"""Unified features for individual sequences, sequence sets, and paired chains."""


import numpy as np
import pandas as pd

from .embeddings import b100
from .utility import frequency_select_physical
from .preprocessing import validate_preprocessing
from ._sequence import _sequences, _length, _attributes, embed_vectors
from .transforms import TransformResult, _transform_vectors, _export_features


class Sequence:
    """Encode a sequence using numeric channels and a selectable transform.

    Parameters
    ----------
    sequence : str
        Nonempty, case-sensitive residue string. All symbols must be present
        in every supplied embedding channel; no unknown-residue imputation is done.
    attributes : mapping, optional
        Metadata exposed as attributes. Names colliding with API or computed
        state are rejected.
    """

    def __init__(self, sequence, attributes=None):
        self.sequence = _sequences([sequence])[0]
        reserved = {"sequence", "embeddings", "frequency", "spectrum_length",
                    "spectra", "cross_spectra", "spectrum_normalization", "spectrum_centered", "result"}
        if attributes is not None:
            for name, value in attributes.items():
                if not isinstance(name, str) or name.startswith('_') or name in reserved or hasattr(self, name):
                    raise ValueError(f"reserved or invalid attribute name: {name!r}")
                setattr(self, name, value)

    def embed_sequence(self, embeddings):
        """Store one numeric vector per channel from a mapping of mappings.

        ``embeddings[channel][residue]`` must be a finite real number for
        every residue in the sequence. Missing symbols raise ValueError.
        Re-embedding invalidates previously computed spectra.
        """
        vectors = embed_vectors(self.sequence, embeddings)
        self.embeddings = vectors
        for name in ("spectra", "cross_spectra", "frequency", "spectrum_length",
                     "spectrum_normalization", "spectrum_centered", "result"):
            self.__dict__.pop(name, None)

    def compute_spectra(self, length, *, normalization="none", center=False):
        """Store frequency-centered FFTs after optional signal scaling and padding.

        Parameters
        ----------
        length : int
            FFT length, at least the sequence length. Each stored spectrum has
            ``length + 1`` entries: the first centered bin is repeated at the
            end for periodic closure. ``frequency`` contains matching cycles
            per residue, including the repeated bin shifted by one cycle.
        normalization : {"none", "length", "sqrt_length", "energy"}, optional
            Divide channel values by 1, residue count N, sqrt(N), or L2 norm,
            respectively, before padding. Energy uses the optionally centered
            channel; zero-energy channels stay zero. Defaults to raw scaling.
        center : bool, optional
            Subtract the channel mean over actual residues before scaling and
            padding. This is distinct from centering FFT bins with fftshift.
        """
        self.compute_transform(length, normalization=normalization, center=center)
        self._store_fourier_compatibility()

    def compute_transform(self, length=None, *, method="fourier", normalization="none",
                          center=False, **options):
        """Store and return a TransformResult from any supported transform.

        Embed first. All backends share centering/normalization, coordinate-aware
        feature export, within-sequence products, and explicit fixed-grid rules.
        Recomputing invalidates Fourier compatibility attributes, even on failure.
        """
        for name in ("result", "spectra", "cross_spectra", "frequency", "spectrum_length",
                     "spectrum_normalization", "spectrum_centered"):
            self.__dict__.pop(name, None)
        if not hasattr(self, "embeddings"):
            raise ValueError("call embed_sequence before compute_transform")
        self.result = _transform_vectors(self.embeddings, len(self.sequence),
            transform=method, length=length, normalization=normalization, center=center, **options)
        return self.result

    def features(self, component="amplitude", *, kinds=("coefficients",), selector=None):
        """Return values/names using the same export interface for every backend."""
        if not hasattr(self, 'result'):
            raise ValueError("call compute_transform before features")
        return self.result.features(component, kinds=kinds, selector=selector)

    def _store_fourier_compatibility(self):
        """Adapt a Fourier result to the historical closed-bin attributes."""
        closed = self.result._with_periodic_closure()
        self.frequency, self.spectra = closed.coordinates['frequency'], closed.coefficients
        self.spectrum_length = self.result.metadata['length']
        self.spectrum_normalization = self.result.metadata['normalization']
        self.spectrum_centered = self.result.metadata['center']

    def _legacy_result(self):
        """Read current compatibility arrays, including explicit caller edits."""
        return TransformResult(self.spectra, {'frequency': self.frequency},
            {'transform': 'fourier', 'length': self.spectrum_length})

    def compute_cross_spectra(self):
        """Store upper-triangle channel products ``X_i * conjugate(X_j)``.

        Includes auto-spectra; keys are ``channel_i_channel_j``. Products
        inherit both input channels' preprocessing. There is no additional
        normalization or statistical averaging of cross-products.
        """
        if not hasattr(self, "spectra"):
            raise ValueError("call compute_spectra before compute_cross_spectra")
        self.cross_spectra = self._legacy_result().cross().coefficients

    def compute_spectral_properties(self, spectra_types=("spectra",),
                                    select_method=frequency_select_physical,
                                    compute_amplitude=True, compute_phase=False):
        """Return selected spectral amplitudes and/or phases as a dictionary.

        Parameters
        ----------
        spectra_types : iterable of str
            ``spectra`` and/or ``cross_spectra``, computed beforehand.
        select_method : callable
            Receives ``frequency`` and returns a same-shape boolean mask.
        compute_amplitude, compute_phase : bool
            Include absolute magnitudes and/or phases in radians, respectively.

        Returns
        -------
        dict
            Scalar features named ``Amp_<channel>_<frequency>`` or
            ``Phase_<channel>_<frequency>``. Magnitudes inherit preprocessing
            from compute_spectra; this export applies no further normalization.
        """
        if not hasattr(self, "spectra"):
            raise ValueError("call compute_spectra before exporting properties")
        datasets = {}
        for kind in spectra_types:
            if kind not in ("spectra", "cross_spectra") or not hasattr(self, kind):
                raise ValueError(f"unavailable spectrum type: {kind!r}")
            if kind in datasets:
                raise ValueError(f'duplicate spectrum type: {kind}')
            datasets[kind] = getattr(self, kind)
        components = tuple(name for name, enabled in
                           [('amplitude', compute_amplitude), ('phase', compute_phase)] if enabled)
        values, names = _export_features(datasets, {'frequency': self.frequency},
            'fourier', components, lambda c: select_method(c['frequency']), legacy=True)
        return dict(zip(names, values))


class SequenceFFT(Sequence):
    """Compatibility wrapper for the historical Fourier sequence API.

    New code can use Sequence.compute_transform for any backend. compute_spectra
    and compute_spectral_properties preserve the original closed bins and names.
    """


def _feature_table(items, component, kinds, selector, attributes):
    """Export generic features with shared coordinates and protected metadata."""
    rows, schema = [], None
    for item in items:
        values, names = item.features(component, kinds=kinds, selector=selector)
        if schema is not None and not np.array_equal(names, schema):
            raise ValueError('selector produced inconsistent columns across sequences')
        schema = names
        row = {name: getattr(item, name) for name in (() if attributes is None else attributes)}
        if row.keys() & set(names):
            raise ValueError('metadata names must not collide with features')
        row.update(zip(names, values))
        rows.append(row)
    return pd.DataFrame(rows)


class SequenceSet:
    """Compute features for sequences sharing one embedding and transform grid.

    Parameters
    ----------
    embedding : mapping, optional
        Channel-to-residue score mappings; defaults to canonical BLOSUM100.
    compute_cross_spectra : bool, optional
        Compute within-sequence channel cross-products, including auto-spectra.
    method : {'fourier', 'z', 'mellin', 'wavelet'}, optional
        Transform used by all sequences. Defaults to Fourier.
    transform_options : mapping, optional
        Backend options such as radius, mellin_c, wavelet, level, and mode.
    normalization : {"none", "length", "sqrt_length", "energy"}, optional
        Per-channel scaling before transformation.
    center : bool, optional
        Center each channel over actual residues before scaling and padding.
    """

    def __init__(self, embedding=b100, compute_cross_spectra=True, *,
                 normalization="none", center=False, method="fourier", transform_options=None):
        self.method, self.transform_options = method, transform_options
        self.embedding = embedding
        self.compute_cross_spectra = compute_cross_spectra
        validate_preprocessing(normalization, center)
        self.normalization = normalization
        self.center = center

    def from_sequence_list(self, sequences, spectrum_length=None, attributes=None, *, length=None):
        """Populate the set from an iterable of nonempty sequence strings.

        ``length`` defaults to the smallest power of two covering the longest
        sequence. ``spectrum_length`` is the compatibility spelling; do not
        provide both. ``attributes`` maps metadata names to equal-length
        columns, one value per sequence. Mutates this object and returns None.
        """
        for name in ('sequences', 'spectrum_length', 'length', 'monomers', 'results'):
            self.__dict__.pop(name, None)
        if length is not None and spectrum_length is not None:
            raise ValueError('use length or the compatibility spectrum_length, not both')
        sequences = _sequences(sequences)
        length = _length(length if length is not None else spectrum_length, max(map(len, sequences)))
        rows = _attributes(attributes, len(sequences))
        monomers = []
        for sequence, row in zip(sequences, rows):
            item_class = SequenceFFT if self.method == 'fourier' else Sequence
            item = item_class(sequence, row)
            item.embed_sequence(self.embedding)
            item.compute_transform(length, method=self.method, normalization=self.normalization,
                                   center=self.center, **(self.transform_options or {}))
            if self.method == 'fourier':
                item._store_fourier_compatibility()
                if self.compute_cross_spectra:
                    item.compute_cross_spectra()
            monomers.append(item)
        self.sequences, self.spectrum_length, self.monomers = sequences, length, monomers
        self.results = [item.result for item in monomers]
        self.length = length

    def export_features(self, component="amplitude", *, kinds=("coefficients",),
                        selector=None, attributes=None):
        """Return a DataFrame for any backend using TransformResult.features.

        All sequences share one grid. selector receives a coordinate dictionary;
        kinds can include aligned within-sequence cross_coefficients. Metadata
        columns must not collide with feature names. Populate the set first.
        """
        if not hasattr(self, 'monomers'):
            raise ValueError('call from_sequence_list before export_features')
        return _feature_table(self.monomers, component, kinds, selector, attributes)

    def export_datatable(self, spectra=("spectra",), attributes=None,
                         select_method=frequency_select_physical,
                         compute_amplitude=True, compute_phase=False):
        """Return a pandas DataFrame with one row per sequence.

        ``spectra`` chooses computed spectrum types. ``attributes`` lists metadata
        columns to copy; selectors and amplitude/phase flags follow
        :meth:`SequenceFFT.compute_spectral_properties`. Populate the set first.
        """
        if not hasattr(self, "monomers"):
            raise ValueError("call from_sequence_list before export_datatable")
        if not hasattr(self.monomers[0], 'frequency'):
            raise ValueError('legacy frequency export is Fourier-only; use export_features')
        return _datatable(self.monomers, spectra, attributes, select_method,
                          compute_amplitude, compute_phase)


def _datatable(items, spectra, attributes, selector, amplitude, phase):
    """Export items without allowing metadata to overwrite spectral features."""
    rows = []
    for item in items:
        row = {name: getattr(item, name) for name in (() if attributes is None else attributes)}
        features = item.compute_spectral_properties(spectra, selector, amplitude, phase)
        if row.keys() & features.keys():
            raise ValueError("metadata names must not collide with spectral features")
        row.update(features)
        rows.append(row)
    return pd.DataFrame(rows)


class ComplexSet:
    """Compute ordered cross-chain products for paired sequences.

    ``embedding`` defaults to BLOSUM100. The historical
    ``compute_cross_cross_spectra`` option is unimplemented; True raises
    NotImplementedError rather than silently doing nothing. ``normalization``
    and ``center`` follow SequenceFFT.compute_spectra and are applied to each
    chain independently before cross-products. method selects fourier, z,
    mellin, or wavelet; transform_options supplies the backend options. Use
    export_features for every backend. export_datatable is a Fourier-only
    compatibility adapter with the original frequency selectors and names.
    """

    def __init__(self, embedding=b100, compute_cross_cross_spectra=False, *,
                 normalization="none", center=False, method="fourier", transform_options=None):
        self.method, self.transform_options = method, transform_options
        if compute_cross_cross_spectra:
            raise NotImplementedError("cross-cross spectra are not implemented")
        self.embedding = embedding
        self.compute_cross_cross_spectra = compute_cross_cross_spectra
        validate_preprocessing(normalization, center)
        self.normalization = normalization
        self.center = center

    def from_sequence_list(self, A_sequences, B_sequences, spectrum_length=None, attributes=None, *, length=None):
        """Populate equally sized collections of paired A/B sequence strings.

        All chains share ``length`` (or compatibility ``spectrum_length``),
        defaulting to the smallest power of two covering the longest chain. Attribute columns have one value per pair.
        ``A`` and ``B`` store chain results; ``AB`` stores every ordered channel
        product after each chain's preprocessing. Fourier compatibility attributes
        are also available when method='fourier'.
        """
        for name in ('A_sequences', 'B_sequences', 'spectrum_length', 'length', 'A', 'B', 'AB', 'results'):
            self.__dict__.pop(name, None)
        if length is not None and spectrum_length is not None:
            raise ValueError('use length or the compatibility spectrum_length, not both')
        a_sequences, b_sequences = _sequences(A_sequences), _sequences(B_sequences)
        if len(a_sequences) != len(b_sequences):
            raise ValueError("A_sequences and B_sequences must have equal lengths")
        length = _length(length if length is not None else spectrum_length,
                         max(map(len, a_sequences + b_sequences)))
        rows = _attributes(attributes, len(a_sequences))
        a_items, b_items, pairs = [], [], []
        for a, b, row in zip(a_sequences, b_sequences, rows):
            item_class = SequenceFFT if self.method == 'fourier' else Sequence
            first, second = item_class(a, row), item_class(b, row)
            for item in (first, second):
                item.embed_sequence(self.embedding)
                item.compute_transform(length, method=self.method, normalization=self.normalization,
                                       center=self.center, **(self.transform_options or {}))
                if self.method == 'fourier':
                    item._store_fourier_compatibility()
            pair = item_class(a, row)
            pair.sequence = (a, b)
            pair.result = first.result.cross(second.result)
            if self.method == 'fourier':
                pair._store_fourier_compatibility()
            a_items.append(first)
            b_items.append(second)
            pairs.append(pair)
        self.A_sequences, self.B_sequences = a_sequences, b_sequences
        self.spectrum_length = length
        self.A, self.B, self.AB = a_items, b_items, pairs
        self.results = [item.result for item in pairs]
        self.length = length

    def export_features(self, component="amplitude", *, chain_ID="AB", kinds=("coefficients",),
                        selector=None, attributes=None):
        """Export either chain or ordered A/B coefficient products for any backend.

        chain_ID is A, B, or AB. Paired products use every A/B channel combination;
        wavelet products remain band/position features, not Fourier cross-spectra.
        """
        if chain_ID not in ('A', 'B', 'AB'):
            raise ValueError("chain_ID must be 'A', 'B', or 'AB'")
        if not hasattr(self, 'AB'):
            raise ValueError('call from_sequence_list before export_features')
        return _feature_table(getattr(self, chain_ID), component, kinds, selector, attributes)

    def export_datatable(self, chain_ID="AB", attributes=None,
                         select_method=frequency_select_physical,
                         compute_amplitude=True, compute_phase=False):
        """Return a feature DataFrame for ``chain_ID`` A, B, or AB.

        ``attributes`` lists metadata columns; frequency selection and feature
        flags follow :meth:`SequenceFFT.compute_spectral_properties`.
        Populate the paired set before calling this method.
        """
        if chain_ID not in ("A", "B", "AB"):
            raise ValueError("chain_ID must be 'A', 'B', or 'AB'")
        if not hasattr(self, "AB"):
            raise ValueError("call from_sequence_list before export_datatable")
        if not hasattr(self.AB[0], 'frequency'):
            raise ValueError('legacy frequency export is Fourier-only; use export_features')
        return _datatable(getattr(self, chain_ID), ("spectra",), attributes,
                          select_method, compute_amplitude, compute_phase)
