"""Optional scikit-learn integration; install Protspect with the ``ml`` extra."""

from copy import deepcopy
from collections import Counter
from collections.abc import Mapping

import numpy as np
from scipy.fft import fftfreq, fftshift

try:
    from sklearn.base import BaseEstimator, TransformerMixin
    from sklearn.utils.validation import check_is_fitted
except ImportError as error:
    raise ImportError("SpectralTransformer requires scikit-learn; install Protspect[ml]") from error

from . import __version__
from .embeddings import b100
from .features import SequenceFFT, _length, _sequences
from .preprocessing import validate_preprocessing
from .utility import frequency_select_physical


def _copy_embedding(embedding):
    """Validate and snapshot every channel, including symbols unseen in training."""
    if not isinstance(embedding, Mapping) or not embedding:
        raise ValueError("embedding must be a nonempty mapping of channels")
    result = {}
    alphabet = None
    for name, values in embedding.items():
        if not isinstance(name, str) or not name or not isinstance(values, Mapping) or not values:
            raise ValueError("embedding channels require nonempty names and residue mappings")
        if any(not isinstance(aa, str) or len(aa) != 1 for aa in values):
            raise ValueError("embedding symbols must be single-character strings")
        if alphabet is not None and set(values) != alphabet:
            raise ValueError("all embedding channels must cover the same alphabet")
        alphabet = set(values)
        if any(np.ndim(value) != 0 or np.iscomplexobj(value) for value in values.values()):
            raise ValueError("embedding values must be finite real scalars")
        try:
            channel = {aa: float(value) for aa, value in values.items()}
        except (TypeError, ValueError, OverflowError) as error:
            raise ValueError("embedding values must be finite real scalars") from error
        if not np.isfinite(list(channel.values())).all():
            raise ValueError("embedding values must be finite real scalars")
        result[name] = channel
    return result


def _channel_statistics(sequences, embedding, standardize):
    """Compute residue-weighted population moments using only training residues."""
    counts = Counter(aa for sequence in sequences for aa in sequence)
    unknown = set(counts) - set(next(iter(embedding.values())))
    if unknown:
        raise ValueError(f"training sequences contain unknown residues: {sorted(unknown)}")
    residues = list(counts)
    weights = np.array([counts[aa] for aa in residues], dtype=float)
    weights /= weights.sum()
    means, scales, effective = {}, {}, {}
    for name, channel in embedding.items():
        values = np.array([channel[aa] for aa in residues])
        mean, scale = 0.0, 1.0
        if standardize:
            peak = np.max(np.abs(values))
            unit = values / peak if peak else np.zeros_like(values)
            unit_mean = float(weights @ unit)
            mean = float(values[0]) if np.all(values == values[0]) else float(peak * unit_mean)
            deviations = unit - (mean / peak if peak else 0.0)
            scale = float(peak * np.sqrt(weights @ deviations**2))
        if not np.isfinite([mean, scale]).all():
            raise ValueError("channel statistics overflowed; rescale the embedding")
        means[name], scales[name] = mean, scale if scale > 0 else 1.0
        effective[name] = {aa: (value - mean) / scales[name] for aa, value in channel.items()}
        if not np.isfinite(list(effective[name].values())).all():
            raise ValueError("channel standardization produced nonfinite values")
    return means, scales, effective


class SpectralTransformer(TransformerMixin, BaseEstimator):
    """Fit a fixed Fourier feature schema for a one-dimensional sequence iterable.

    Parameters
    ----------
    embedding : dict of dict, optional
        Numeric channels indexed by single-character residue codes. All channels
        must cover the same alphabet. None uses BLOSUM100. For estimator cloning,
        supply ordinary dictionaries for custom embeddings.
    spectrum_length : int, optional
        FFT length. None selects the smallest power of two covering the longest
        training sequence. Transform rejects sequences longer than this fitted
        grid. No cropping or grid expansion is performed.
    normalization : {"none", "length", "sqrt_length", "energy"}, optional
        Per-sequence scaling before padding; see SequenceFFT.compute_spectra.
    center : bool, optional
        Subtract each channel's mean over the current sequence before scaling.
    standardize_channels : bool, optional
        Learn a mean and population standard deviation for each channel across
        all training residues, with each residue weighted equally. Apply these
        frozen statistics before per-sequence centering/scaling. Zero-variance
        channels use scale 1. Defaults to False.
    spectra_types : tuple of str, optional
        Export ``spectra`` and/or ``cross_spectra``; defaults to spectra only.
    select_method : callable, optional
        Boolean frequency selector, evaluated once during fit. The resulting
        mask is reused unchanged, including for a stateful callable.
    compute_amplitude, compute_phase : bool, optional
        Export magnitudes (default True) and/or angles in radians (default False).

    Attributes
    ----------
    spectrum_length_ : int
        Fitted grid length, independent of later batch lengths.
    frequency_, selected_frequencies_ : numpy.ndarray
        Full grid including its periodic closure and the selected grid subset.
    feature_names_out_ : numpy.ndarray of str
        Frozen output column order, also available via get_feature_names_out.
    channel_mean_, channel_scale_ : dict
        Fitted channel offsets/scales, or 0/1 when standardization is disabled.
    embedding_ : dict of dict
        Independent snapshot of the original embedding before standardization.
    n_features_out_ : int
        Number of output features.

    Notes
    -----
    Put this transformer inside a Pipeline before cross-validation so grid and
    channel statistics are fitted only on each training fold. Labels are ignored.
    Output is a float64 array; ``set_output(transform="pandas")`` is also supported.
    Use set_params followed by fit to change configuration. Direct edits to
    constructor parameters take effect on the next fit; transform uses snapshots.
    This transformer handles single sequences, not paired chains.
    """

    def __init__(self, *, embedding=None, spectrum_length=None, normalization="none",
                 center=False, standardize_channels=False, spectra_types=("spectra",),
                 select_method=frequency_select_physical, compute_amplitude=True,
                 compute_phase=False):
        self.embedding = embedding
        self.spectrum_length = spectrum_length
        self.normalization = normalization
        self.center = center
        self.standardize_channels = standardize_channels
        self.spectra_types = spectra_types
        self.select_method = select_method
        self.compute_amplitude = compute_amplitude
        self.compute_phase = compute_phase

    def _clear_fit(self):
        """Invalidate fitted snapshots without altering constructor parameters."""
        for name in list(vars(self)):
            if name.endswith('_') and not name.startswith('_'):
                delattr(self, name)

    def set_params(self, **params):
        """Set estimator parameters and invalidate any previously fitted schema."""
        result = super().set_params(**params)
        if params:
            self._clear_fit()
        return result

    def fit(self, X, y=None):
        """Freeze grid, selector, embedding, and channel moments from X alone.

        X is a nonempty iterable of nonempty sequence strings. y is accepted for
        Pipeline compatibility and ignored. Return this fitted transformer.
        """
        self._clear_fit()
        sequences = _sequences(X)
        validate_preprocessing(self.normalization, self.center)
        for name in ('standardize_channels', 'compute_amplitude', 'compute_phase'):
            if not isinstance(getattr(self, name), (bool, np.bool_)):
                raise ValueError(f"{name} must be a boolean")
        if not self.compute_amplitude and not self.compute_phase:
            raise ValueError("enable amplitude and/or phase features")
        if not isinstance(self.spectra_types, (list, tuple)) or not self.spectra_types:
            raise ValueError("spectra_types must be a nonempty list or tuple")
        kinds = tuple(self.spectra_types)
        if any(kind not in ('spectra', 'cross_spectra') for kind in kinds) or len(set(kinds)) != len(kinds):
            raise ValueError("spectra_types must uniquely select spectra and/or cross_spectra")
        length = _length(self.spectrum_length, max(map(len, sequences)))
        embedding = _copy_embedding(b100 if self.embedding is None else self.embedding)
        means, scales, effective = _channel_statistics(sequences, embedding, self.standardize_channels)
        bins = fftshift(fftfreq(length))
        frequencies = np.append(bins, bins[0] + 1)
        if not callable(self.select_method):
            raise ValueError("select_method must be callable")
        mask = np.asarray(self.select_method(frequencies.copy()))
        if mask.dtype != bool or mask.shape != frequencies.shape or not mask.any():
            raise ValueError("select_method must select at least one bin with a same-shape boolean mask")
        # Build fitted state on a temporary object, so failed fits leave no usable
        # old or partially updated model behind.
        candidate = type(self)()
        candidate.package_version_ = __version__
        candidate.standardize_channels_ = bool(self.standardize_channels)
        candidate.spectrum_length_ = length
        candidate.frequency_ = frequencies.copy()
        candidate.selected_mask_ = mask.copy()
        candidate.selected_frequencies_ = frequencies[mask].copy()
        candidate.embedding_ = embedding
        candidate.effective_embedding_ = effective
        candidate.channel_mean_, candidate.channel_scale_ = means, scales
        candidate.normalization_, candidate.center_ = self.normalization, bool(self.center)
        candidate.spectra_types_ = kinds
        candidate.compute_amplitude_ = bool(self.compute_amplitude)
        candidate.compute_phase_ = bool(self.compute_phase)
        first = candidate._row(sequences[0])
        candidate.feature_names_out_ = np.asarray(list(first), dtype=object)
        candidate.n_features_out_ = len(first)
        self.__dict__.update({name: value for name, value in vars(candidate).items()
                              if name.endswith('_') and not name.startswith('_')})
        return self

    def _row(self, sequence):
        """Compute one feature dictionary using only frozen fitted settings."""
        item = SequenceFFT(sequence)
        item.embed_sequence(self.effective_embedding_)
        item.compute_spectra(self.spectrum_length_, normalization=self.normalization_, center=self.center_)
        if 'cross_spectra' in self.spectra_types_:
            item.compute_cross_spectra()
        return item.compute_spectral_properties(self.spectra_types_,
            lambda frequencies: self.selected_mask_, self.compute_amplitude_, self.compute_phase_)

    def transform(self, X):
        """Return fixed-width features without updating fitted statistics or bins.

        Unknown residues, empty sequences, and sequences exceeding the fitted
        grid raise ValueError. transform never refits or silently truncates.
        """
        check_is_fitted(self, ['feature_names_out_', 'spectrum_length_'])
        sequences = _sequences(X)
        if max(map(len, sequences)) > self.spectrum_length_:
            raise ValueError("sequence exceeds fitted spectrum_length_; choose a larger explicit grid and refit")
        rows = []
        for sequence in sequences:
            row = self._row(sequence)
            if list(row) != self.feature_names_out_.tolist():
                raise ValueError("computed feature columns do not match the fitted schema")
            rows.append(list(row.values()))
        values = np.asarray(rows, dtype=float)
        if not np.isfinite(values).all():
            raise ValueError("spectral features must be finite")
        return values

    def fit_transform(self, X, y=None, **fit_params):
        """Fit and transform, materializing generators only once."""
        if fit_params:
            raise TypeError("SpectralTransformer accepts no additional fit parameters")
        sequences = _sequences(X)
        return self.fit(sequences, y).transform(sequences)

    def get_feature_names_out(self, input_features=None):
        """Return a copy of the fitted output names in transform column order."""
        check_is_fitted(self, ['feature_names_out_'])
        if input_features is not None and len(input_features) != 1:
            raise ValueError("SpectralTransformer accepts one sequence input per sample")
        return self.feature_names_out_.copy()

    def get_feature_metadata(self):
        """Return independent, JSON-serializable fitted feature provenance.

        Includes exact embedding scores, fitted channel moments, full and
        selected grids, the selector mask, normalization options, and column
        order. Store this beside exported features. It describes the fitted
        representation; use pickle/joblib to persist the executable Pipeline.
        """
        check_is_fitted(self, ['feature_names_out_'])
        return {
            "package": "Protspect", "package_version": self.package_version_,
            "spectrum_length": self.spectrum_length_,
            "normalization": self.normalization_, "center": self.center_,
            "standardize_channels": self.standardize_channels_,
            "channel_weighting": "one weight per training residue; padding excluded",
            "channel_mean": dict(self.channel_mean_), "channel_scale": dict(self.channel_scale_),
            "embedding": deepcopy(self.embedding_),
            "frequency": self.frequency_.tolist(),
            "selected_mask": self.selected_mask_.tolist(),
            "selected_frequencies": self.selected_frequencies_.tolist(),
            "spectra_types": list(self.spectra_types_),
            "compute_amplitude": self.compute_amplitude_, "compute_phase": self.compute_phase_,
            "feature_names": self.get_feature_names_out().tolist(),
        }


class SequenceTransformer(TransformerMixin, BaseEstimator):
    """Fixed-width Fourier, Z, Mellin, or wavelet features for sklearn pipelines.

    ``method`` selects fourier, z, mellin, or wavelet; ``component`` selects
    amplitude, phase, real, or imaginary output. ``embedding`` defaults to BLOSUM100.
    ``normalization`` and ``center`` act on actual residues before padding.
    ``transform_options`` is a mapping of keyword arguments to
    :func:`protspect.transform_sequence` (radius, mellin_c, wavelet, level, mode).
    ``length`` is learned from training sequences alone when omitted. Longer
    held-out sequences raise rather than changing the schema. Channel moments
    are also learned from training residues alone when standardization is enabled.
    Constructor mutations take effect only on refit; fitted settings are copied.
    """

    def __init__(self, *, embedding=None, method="fourier", length=None,
                 normalization="none", center=False, standardize_channels=False,
                 component="amplitude", transform_options=None):
        self.embedding = embedding
        self.method = method
        self.length = length
        self.normalization = normalization
        self.center = center
        self.standardize_channels = standardize_channels
        self.component = component
        self.transform_options = transform_options

    def set_params(self, **params):
        """Invalidate fitted settings when estimator parameters change."""
        result = super().set_params(**params)
        if params:
            self._clear_fit()
        return result

    def _clear_fit(self):
        """Remove fitted state, including after an unsuccessful refit."""
        for name in list(vars(self)):
            if name.endswith('_') and not name.startswith('_'):
                delattr(self, name)

    def fit(self, X, y=None):
        """Freeze training-only length, embedding, options, and output schema."""
        from .transforms import transform_sequence

        self._clear_fit()
        sequences = _sequences(X)
        length = _length(self.length, max(map(len, sequences)))
        if not isinstance(self.standardize_channels, (bool, np.bool_)):
            raise ValueError("standardize_channels must be boolean")
        embedding = _copy_embedding(b100 if self.embedding is None else self.embedding)
        means, scales, effective = _channel_statistics(sequences, embedding, self.standardize_channels)
        options = {} if self.transform_options is None else deepcopy(self.transform_options)
        allowed = {'radius', 'mellin_c', 'wavelet', 'level', 'mode'}
        if not isinstance(options, dict) or not set(options) <= allowed:
            raise ValueError(f"transform_options must be a mapping using {sorted(allowed)}")
        settings = dict(transform=self.method, length=length,
                        normalization=self.normalization, center=self.center, **options)
        names = None
        for sequence in sequences:
            result = transform_sequence(sequence, effective, **settings)
            _, current_names = result.features(self.component)
            if names is not None and not np.array_equal(names, current_names):
                raise ValueError("training sequences produced inconsistent feature schemas")
            names = current_names
        self.length_ = length
        self.embedding_, self.effective_embedding_ = embedding, effective
        self.channel_mean_, self.channel_scale_ = means, scales
        self.settings_, self.component_ = deepcopy(settings), self.component
        self.feature_names_out_ = names.copy()
        self.n_features_out_ = len(names)
        # Metadata excludes training-sequence-specific residue count.
        self.transform_metadata_ = {k: deepcopy(v) for k, v in result.metadata.items()
                                    if k != 'sequence_length'}
        self.coordinates_ = deepcopy(result.coordinates)
        return self

    def transform(self, X):
        """Apply only frozen fitted settings; return channel-major float features."""
        from .transforms import transform_sequence

        check_is_fitted(self, ['feature_names_out_', 'settings_'])
        sequences = _sequences(X)
        rows = []
        for sequence in sequences:
            result = transform_sequence(sequence, self.effective_embedding_, **self.settings_)
            values, names = result.features(self.component_)
            if not np.array_equal(names, self.feature_names_out_):
                raise ValueError("computed features do not match the fitted schema")
            rows.append(values)
        return np.asarray(rows, dtype=float)

    def fit_transform(self, X, y=None, **fit_params):
        """Materialize generators once before fitting and transforming."""
        if fit_params:
            raise TypeError("fit parameters are not supported")
        sequences = _sequences(X)
        return self.fit(sequences, y).transform(sequences)

    def get_feature_names_out(self, input_features=None):
        """Return fitted feature names in output column order."""
        check_is_fitted(self, ['feature_names_out_'])
        return self.feature_names_out_.copy()
