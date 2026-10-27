"""Optional scikit-learn integration; install Protspect with the ``ml`` extra."""

from copy import deepcopy
from collections import Counter
from collections.abc import Mapping

import numpy as np

try:
    from sklearn.base import BaseEstimator, TransformerMixin
    from sklearn.utils.validation import check_is_fitted
except ImportError as error:
    raise ImportError("Protspect sklearn integration requires scikit-learn; install Protspect[ml]") from error

from . import __version__
from .embeddings import b100
from ._sequence import _length, _sequences
from .transforms import transform_sequence, transform_pair, _selection, _export_features
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


def _samples(X, paired):
    """Materialize single sequences or exact two-chain samples without truncation."""
    if not paired:
        return _sequences(X)
    if isinstance(X, str):
        raise TypeError('provide an iterable of two-sequence pairs')
    samples = []
    for pair in X:
        if isinstance(pair, str):
            raise ValueError('each paired sample must contain exactly two sequences')
        pair = list(pair)
        if len(pair) != 2:
            raise ValueError('each paired sample must contain exactly two sequences')
        samples.append(tuple(_sequences(pair)))
    if not samples:
        raise ValueError('provide at least one two-sequence pair')
    return samples


def _jsonable(value):
    """Copy numerical provenance into JSON-safe values, including Z-plane axes."""
    if isinstance(value, np.ndarray):
        if np.iscomplexobj(value):
            return {'real': value.real.tolist(), 'imaginary': value.imag.tolist()}
        return value.tolist()
    if isinstance(value, np.generic):
        return _jsonable(value.item())
    if isinstance(value, complex):
        return {'real': value.real, 'imaginary': value.imag}
    if isinstance(value, dict):
        return {k: _jsonable(v) for k,v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(v) for v in value]
    return deepcopy(value)


class SequenceTransformer(TransformerMixin, BaseEstimator):
    """Fit one coefficient-feature schema for every transform and single/paired inputs.

    method selects fourier, z, mellin, or wavelet. embedding defaults to BLOSUM100.
    length defaults to the next power of two covering training chains only.
    Longer held-out chains raise. normalization and center operate on each real
    chain before padding. standardize_channels learns training-residue moments.

    component is amplitude, phase, real, imaginary, or an ordered tuple thereof.
    kinds selects coefficients and/or cross_coefficients. selector receives a
    dictionary of coordinate arrays and returns a same-shape boolean mask; an
    explicit boolean mask is also accepted. It is evaluated only at fit, copied,
    and frozen. Transform-specific options (radius, mellin_c, wavelet, level,
    mode) go in transform_options. With paired=True, X contains (A, B) samples
    and features describe all ordered A/B channel products. Non-Fourier products
    remain coefficient products rather than frequency-domain cross-spectra.

    Use inside a Pipeline so grid and channel moments fit within each CV fold.
    get_feature_metadata returns independent JSON-safe provenance. set_output
    supports pandas. Constructor mutations apply only after refitting; set_params
    invalidates fit. The SpectralTransformer compatibility subclass delegates
    here while preserving the original Fourier names, selection, and closure.
    """

    def __init__(self, *, embedding=None, method="fourier", length=None,
                 normalization="none", center=False, standardize_channels=False,
                 component="amplitude", transform_options=None, selector=None,
                 kinds=("coefficients",), paired=False):
        self.embedding = embedding
        self.method = method
        self.length = length
        self.normalization = normalization
        self.center = center
        self.standardize_channels = standardize_channels
        self.component = component
        self.transform_options = transform_options
        self.selector = selector
        self.kinds = kinds
        self.paired = paired

    def __setstate__(self, state):
        """Upgrade fitted snapshots saved by the preceding two-interface API."""
        super().__setstate__(state)
        if 'method' in vars(self):
            for name, value in [('selector', None), ('kinds', ('coefficients',)), ('paired', False)]:
                if not hasattr(self, name):
                    setattr(self, name, value)
        if not hasattr(self, 'feature_names_out_'):
            return
        if not hasattr(self, 'settings_') and hasattr(self, 'spectrum_length_'):
            self.length_ = self.spectrum_length_
            self.settings_ = dict(method='fourier', length=self.length_,
                                  normalization=self.normalization_, center=self.center_)
            self.component_ = tuple(c for c, enabled in
                [('amplitude', self.compute_amplitude_), ('phase', self.compute_phase_)] if enabled)
            self.kinds_ = tuple('coefficients' if k == 'spectra' else 'cross_coefficients'
                                for k in self.spectra_types_)
            self.coordinates_ = {'frequency': self.frequency_.copy()}
            self.legacy_, self.paired_ = True, False
            self.transform_metadata_ = dict(transform='fourier', length=self.length_,
                normalization=self.normalization_, center=self.center_)
        elif hasattr(self, 'settings_') and not hasattr(self, 'selected_mask_'):
            count = len(next(iter(self.coordinates_.values())))
            self.selected_mask_ = np.ones(count, dtype=bool)
            self.kinds_, self.legacy_, self.paired_ = ('coefficients',), False, False
            self.standardize_channels_ = bool(self.standardize_channels)
            self.normalization_ = self.settings_['normalization']
            self.center_ = self.settings_['center']
            self.package_version_ = __version__

    def _configuration(self):
        """Translate constructor settings; compatibility subclasses only adapt this."""
        return dict(method=self.method, length=self.length, normalization=self.normalization,
            center=self.center, standardize=self.standardize_channels, component=self.component,
            options=self.transform_options, selector=self.selector, kinds=self.kinds,
            paired=self.paired, legacy=False, embedding=self.embedding)

    def set_params(self, **params):
        """Change settings and invalidate the fitted schema until the next fit."""
        result = super().set_params(**params)
        if params:
            self._clear_fit()
        return result

    def _clear_fit(self):
        """Remove old fitted state before every fit, including unsuccessful refits."""
        for name in list(vars(self)):
            if name.endswith('_') and not name.startswith('_'):
                delattr(self, name)

    @staticmethod
    def _result(sample, embedding, settings, paired, legacy):
        """Use the common numerical engine, with an optional old closed-bin view."""
        result = (transform_pair(*sample, embedding, **settings) if paired
                  else transform_sequence(sample, embedding, **settings))
        return result._with_periodic_closure() if legacy else result

    @staticmethod
    def _features(result, component, kinds, selector, legacy):
        """Use the shared export engine for both primary and compatibility names."""
        if not legacy:
            return result.features(component, kinds=kinds, selector=selector)
        datasets = {kind: (result.coefficients if kind == 'coefficients'
                           else result.cross().coefficients) for kind in kinds}
        return _export_features(datasets, result.coordinates, 'fourier', component,
                                selector, legacy=True)

    def fit(self, X, y=None):
        """Freeze training-only grid, channel moments, selection, and column order."""
        self._clear_fit()
        config = self._configuration()
        validate_preprocessing(config['normalization'], config['center'])
        if any(not isinstance(config[k], (bool, np.bool_)) for k in ('standardize', 'paired')):
            raise ValueError('standardize_channels and paired must be booleans')
        samples = _samples(X, config['paired'])
        chains = [s for pair in samples for s in pair] if config['paired'] else samples
        length = _length(config['length'], max(map(len, chains)))
        embedding = _copy_embedding(b100 if config['embedding'] is None else config['embedding'])
        means, scales, effective = _channel_statistics(chains, embedding, config['standardize'])
        options = {} if config['options'] is None else deepcopy(config['options'])
        allowed = {'radius', 'mellin_c', 'wavelet', 'level', 'mode'}
        if not isinstance(options, Mapping) or not set(options) <= allowed:
            raise ValueError(f'transform_options must be a mapping using {sorted(allowed)}')
        settings = dict(method=config['method'], length=length,
            normalization=config['normalization'], center=bool(config['center']), **options)
        result = self._result(samples[0], effective, settings, config['paired'], config['legacy'])
        mask = _selection(result.coordinates, config['selector'])
        if not mask.any():
            raise ValueError('selector must select at least one coefficient')
        component, kinds = deepcopy(config['component']), deepcopy(config['kinds'])
        values, names = self._features(result, component, kinds, mask, config['legacy'])
        if not names.size or not np.isfinite(values).all():
            raise ValueError('select at least one finite output feature')
        # Commit only after every configuration/schema check has succeeded.
        self.length_, self.settings_ = length, deepcopy(settings)
        self.embedding_, self.effective_embedding_ = embedding, effective
        self.channel_mean_, self.channel_scale_ = means, scales
        self.component_, self.kinds_ = component, kinds
        self.selected_mask_ = mask.copy()
        self.paired_, self.legacy_ = bool(config['paired']), config['legacy']
        self.coordinates_ = deepcopy(result.coordinates)
        self.transform_metadata_ = {k: deepcopy(v) for k,v in result.metadata.items()
                                    if k not in ('sequence_length', 'sequence_lengths')}
        self.feature_names_out_, self.n_features_out_ = names.copy(), len(names)
        self.package_version_ = __version__
        self.standardize_channels_ = bool(config['standardize'])
        self.normalization_, self.center_ = config['normalization'], bool(config['center'])
        if config['legacy']:
            self.spectrum_length_ = length
            self.frequency_ = self.coordinates_['frequency'].copy()
            self.selected_frequencies_ = self.frequency_[mask].copy()
            self.spectra_types_ = tuple('spectra' if k == 'coefficients' else 'cross_spectra'
                                       for k in kinds)
            self.compute_amplitude_ = 'amplitude' in component
            self.compute_phase_ = 'phase' in component
        return self

    def transform(self, X):
        """Apply only frozen settings without changing grids, moments, or selection."""
        check_is_fitted(self, ['feature_names_out_', 'settings_'])
        samples = _samples(X, self.paired_)
        chains = [s for pair in samples for s in pair] if self.paired_ else samples
        if max(map(len, chains)) > self.length_:
            raise ValueError('sequence exceeds fitted length; choose a larger explicit grid and refit')
        rows = []
        for sample in samples:
            result = self._result(sample, self.effective_embedding_, self.settings_,
                                  self.paired_, self.legacy_)
            values, names = self._features(result, self.component_, self.kinds_,
                                           self.selected_mask_, self.legacy_)
            if not np.array_equal(names, self.feature_names_out_):
                raise ValueError('computed features do not match fitted schema')
            rows.append(values)
        return np.asarray(rows, dtype=float)

    def fit_transform(self, X, y=None, **fit_params):
        """Materialize generators only once, then use the same fit/transform path."""
        if fit_params:
            raise TypeError('no additional fit parameters are supported')
        self._clear_fit()
        samples = _samples(X, self._configuration()['paired'])
        return self.fit(samples, y).transform(samples)

    def get_feature_names_out(self, input_features=None):
        """Return independent feature names in output-column order."""
        check_is_fitted(self, ['feature_names_out_'])
        if input_features is not None and len(input_features) != (2 if self.paired_ else 1):
            raise ValueError('input_features must match the single/paired sequence input')
        return self.feature_names_out_.copy()

    def get_feature_metadata(self):
        """Return independent JSON-safe grids, scores, moments, controls, and names."""
        check_is_fitted(self, ['feature_names_out_'])
        return _jsonable(dict(package='Protspect', package_version=self.package_version_,
            settings=self.settings_, paired=self.paired_, component=self.component_, kinds=self.kinds_,
            standardize_channels=self.standardize_channels_, embedding=self.embedding_,
            channel_weighting='one weight per training residue; padding excluded',
            channel_mean=self.channel_mean_, channel_scale=self.channel_scale_,
            transform_metadata=self.transform_metadata_, coordinates=self.coordinates_,
            selected_mask=self.selected_mask_, feature_names=self.feature_names_out_))


class SpectralTransformer(SequenceTransformer):
    """Compatibility adapter for the original fitted Fourier feature API.

    All fitting, transformation, channel preprocessing, and feature extraction
    delegate to SequenceTransformer. The original constructor, closed frequency
    grid, ``Amp_``/``Phase_`` feature names, and get_feature_metadata schema are preserved.
    New code should use SequenceTransformer for every transform.
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

    def _configuration(self):
        """Translate old argument names and frequency selectors into the shared API."""
        for name in ('compute_amplitude', 'compute_phase'):
            if not isinstance(getattr(self, name), (bool, np.bool_)):
                raise ValueError(f'{name} must be a boolean')
        if not self.compute_amplitude and not self.compute_phase:
            raise ValueError('enable amplitude and/or phase features')
        kinds = self.spectra_types
        if (not isinstance(kinds, (list, tuple)) or not kinds or
                len(set(kinds)) != len(kinds) or
                any(k not in ('spectra', 'cross_spectra') for k in kinds)):
            raise ValueError('spectra_types must uniquely select spectra and/or cross_spectra')
        if not callable(self.select_method):
            raise ValueError('select_method must be callable')
        return dict(method='fourier', length=self.spectrum_length,
            normalization=self.normalization, center=self.center,
            standardize=self.standardize_channels, embedding=self.embedding,
            component=tuple(c for c, enabled in [('amplitude', self.compute_amplitude),
                ('phase', self.compute_phase)] if enabled), options=None,
            selector=lambda coordinates: self.select_method(coordinates['frequency']),
            kinds=tuple('coefficients' if k == 'spectra' else 'cross_coefficients' for k in kinds),
            paired=False, legacy=True)

    def get_feature_metadata(self):
        """Preserve the original independent JSON-serializable Fourier provenance."""
        check_is_fitted(self, ['feature_names_out_'])
        return {
            'package': 'Protspect', 'package_version': self.package_version_,
            'spectrum_length': self.spectrum_length_,
            'normalization': self.normalization_, 'center': self.center_,
            'standardize_channels': self.standardize_channels_,
            'channel_weighting': 'one weight per training residue; padding excluded',
            'channel_mean': dict(self.channel_mean_), 'channel_scale': dict(self.channel_scale_),
            'embedding': deepcopy(self.embedding_), 'frequency': self.frequency_.tolist(),
            'selected_mask': self.selected_mask_.tolist(),
            'selected_frequencies': self.selected_frequencies_.tolist(),
            'spectra_types': list(self.spectra_types_),
            'compute_amplitude': self.compute_amplitude_, 'compute_phase': self.compute_phase_,
            'feature_names': self.get_feature_names_out().tolist(),
        }
