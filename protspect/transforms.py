"""Alternative transforms with explicit coefficient coordinates and conventions."""

from dataclasses import dataclass
from copy import deepcopy
from numbers import Integral

import numpy as np
from scipy.fft import fft, fftfreq, fftshift

from .embeddings import b100
from ._sequence import _length, _sequences, embed_vectors
from .preprocessing import prepare_signal, validate_preprocessing


@dataclass
class TransformResult:
    """Channel coefficients and shared axes; no periodic-closure duplicate bin.

    ``coordinates`` maps axis names to coefficient-length arrays. Wavelet band
    and within-band index are coefficient coordinates, not frequencies. Metadata
    records preprocessing and transform conventions. Arrays belong to this result.
    """

    coefficients: dict
    coordinates: dict
    metadata: dict

    def features(self, component="amplitude", *, kinds=("coefficients",), selector=None):
        """Export channel-major values/names for any transform.

        component is amplitude, phase, real, imaginary, or an ordered tuple of
        these. kinds selects coefficients and/or cross_coefficients (upper
        triangle, including auto-products). selector is a boolean mask or a
        callable receiving independent coordinate arrays in a dictionary. For
        example, lambda c: c['frequency'] >= 0 for Fourier/Z, or
        lambda c: c['band'] == 'D1' for wavelets. No selector means all bins.
        Cross-products outside Fourier are aligned coefficient products, not
        estimates of frequency-domain cross-spectral density.
        """
        if isinstance(component, str):
            component = (component,)
        if not isinstance(component, (list, tuple)) or not component:
            raise ValueError("component must select at least one component")
        if not isinstance(kinds, (list, tuple)) or not kinds or len(set(kinds)) != len(kinds):
            raise ValueError("kinds must uniquely select coefficients and/or cross_coefficients")
        datasets = {}
        for kind in kinds:
            if kind == 'coefficients':
                datasets[kind] = self.coefficients
            elif kind == 'cross_coefficients':
                datasets[kind] = self.cross().coefficients
            else:
                raise ValueError("kinds must select coefficients and/or cross_coefficients")
        return _export_features(datasets, self.coordinates, self.metadata['transform'],
                                component, selector)

    def cross(self, other=None):
        """Return aligned conjugate products with explicit within/paired provenance.

        With no other result, take the upper triangle including auto-products.
        With another result, take every ordered A/B channel product. Transform
        settings and coordinate grids must match; residue counts may differ.
        Products are not averaged. Their interpretation depends on the backend.
        """
        paired = other is not None
        other = self if other is None else other
        if not isinstance(other, TransformResult):
            raise TypeError("other must be a TransformResult")
        ignored = {'sequence_length', 'sequence_lengths'}
        first_settings = {k:v for k,v in self.metadata.items() if k not in ignored}
        other_settings = {k:v for k,v in other.metadata.items() if k not in ignored}
        if first_settings != other_settings or self.coordinates.keys() != other.coordinates.keys():
            raise ValueError("cross-products require matching transform settings and coordinates")
        if any(not np.array_equal(v, other.coordinates[k]) for k,v in self.coordinates.items()):
            raise ValueError("cross-products require matching transform settings and coordinates")
        products = {}
        names = list(other.coefficients)
        with np.errstate(over='ignore', invalid='ignore'):
            for i, (first, values) in enumerate(self.coefficients.items()):
                for second in names if paired else names[i:]:
                    key = f'A-{first}_B-{second}' if paired else f'{first}_{second}'
                    if key in products:
                        raise ValueError(f'ambiguous cross-product name: {key}')
                    products[key] = values * np.conjugate(other.coefficients[second])
                    if not np.isfinite(products[key]).all():
                        raise ValueError('cross-product produced nonfinite values; rescale channels')
        metadata = deepcopy(self.metadata)
        metadata['product'] = 'paired' if paired else 'within'
        if paired:
            metadata.pop('sequence_length', None)
            metadata['sequence_lengths'] = [self.metadata.get('sequence_length'),
                                             other.metadata.get('sequence_length')]
        return TransformResult(products, deepcopy(self.coordinates), metadata)

    def _with_periodic_closure(self):
        """Create the legacy Fourier view without recomputing the transform."""
        if self.metadata['transform'] != 'fourier':
            raise ValueError('legacy periodic closure is Fourier-only')
        coordinates = deepcopy(self.coordinates)
        f = coordinates['frequency']
        coordinates['frequency'] = np.append(f, f[0] + 1)
        return TransformResult({k: np.append(v, v[0]) for k,v in self.coefficients.items()},
                               coordinates, deepcopy(self.metadata))


def _selection(coordinates, selector):
    """Resolve and validate a mask without exposing mutable coordinate state."""
    count = len(next(iter(coordinates.values())))
    if any(np.asarray(v).shape != (count,) for v in coordinates.values()):
        raise ValueError('coordinate arrays must be one-dimensional and equal length')
    if selector is None:
        return np.ones(count, dtype=bool)
    selected = selector(deepcopy(coordinates)) if callable(selector) else selector
    selected = np.asarray(selected)
    if selected.dtype != bool or selected.shape != (count,):
        raise ValueError('selector must return a same-shape boolean mask matching coordinates')
    return selected.copy()


def _export_features(datasets, coordinates, method, components, selector=None, *, legacy=False):
    """One feature-export engine for canonical APIs and Fourier adapters."""
    functions = {'amplitude': np.abs, 'phase': np.angle, 'real': np.real, 'imaginary': np.imag}
    if any(c not in functions for c in components) or len(set(components)) != len(components):
        raise ValueError(f'components must uniquely select from {tuple(functions)}')
    selected = _selection(coordinates, selector)
    indices = np.flatnonzero(selected)
    values, names = [], []
    for kind, dataset in datasets.items():
        for channel, coefficients in dataset.items():
            if np.asarray(coefficients).shape != selected.shape or not np.isfinite(coefficients).all():
                raise ValueError('coefficients must be finite and match coordinates')
            for component in components:
                with np.errstate(over='ignore', invalid='ignore'):
                    values.extend(functions[component](coefficients)[selected])
                if legacy:
                    prefix = {'amplitude': 'Amp', 'phase': 'Phase'}[component]
                    names.extend(f'{prefix}_{channel}_{f}' for f in coordinates['frequency'][selected])
                else:
                    part = component if kind == 'coefficients' else 'cross_' + component
                    names.extend(f'{channel}:{method}:{part}:{i}' for i in indices)
    if len(set(names)) != len(names):
        raise ValueError('duplicate feature name; choose unambiguous channel names')
    output = np.asarray(values, dtype=float)
    if not np.isfinite(output).all():
        raise ValueError('feature extraction produced nonfinite values; rescale channels')
    return output, np.asarray(names, dtype=object)


def transform_sequence(sequence, embedding=b100, *, method="fourier", transform=None, length=None,
                       normalization="none", center=False, radius=1.0,
                       mellin_c=0.0, wavelet="db2", level=None, mode="symmetric"):
    """Embed a protein sequence and transform each channel on a fixed grid.

    Parameters
    ----------
    sequence : str
        Nonempty residue sequence; all residues must occur in every channel.
    embedding : mapping
        Channel-to-residue mappings, including AAindex and Kidera mappings.
    method : {'fourier', 'z', 'mellin', 'wavelet'}
        Fourier uses shifted DFT bins. Z evaluates sum x[n] z**(-n) at
        z=radius*exp(2*pi*i*f), n starting at zero. Mellin uses the finite sum
        sum x[n]*(n+1)**(c-1-2*pi*i*nu), n starting at zero. It is a unit-spacing
        quadrature convention on positive positions, not an exact continuous
        transform or an assertion of sequence-length invariance.
    length : int, optional
        Shared padded signal length, at least the residue count; defaults to
        the next power of two. Never truncate sequences. Mellin uses only real
        residue positions; length fixes its common log-frequency grid.
    normalization : str
        'none', 'length', 'sqrt_length', or 'energy', applied before padding.
    center : bool
        Subtract each channel's mean over real residues before normalization.
    radius : float
        Positive finite Z-plane radius. One reproduces Fourier; smaller radii
        emphasize later residues and may overflow on long sequences.
    mellin_c : float
        Real part of the Mellin contour. Positive positions are measured from
        the N terminus in residue units; reversing/shifting changes results.
    wavelet : str
        Discrete PyWavelets wavelet name; requires Protspect[wavelets].
    level : int, optional
        Decomposition depth, bounded by dwt_max_level of the padded grid.
        None uses that maximum; zero returns the undecomposed padded signal.
    mode : str
        PyWavelets boundary extension mode, applied after explicit zero-padding.

    Returns
    -------
    TransformResult
        Raw coefficients, coordinates, and provenance. features exports selected
        coefficients and/or within-channel products. cross(other) forms ordered
        paired products. No selector is implicitly applied. The previous
        ``transform=`` spelling remains a compatibility alias for method.
    """
    transform = _resolve_method(method, transform)
    validate_preprocessing(normalization, center)
    sequence = _sequences([sequence])[0]
    vectors = embed_vectors(sequence, embedding)
    return _transform_vectors(vectors, len(sequence), transform=transform, length=length,
        normalization=normalization, center=center, radius=radius, mellin_c=mellin_c,
        wavelet=wavelet, level=level, mode=mode)


def _transform_vectors(vectors, sequence_length, *, transform="fourier", length=None,
                       normalization="none", center=False, radius=1.0,
                       mellin_c=0.0, wavelet="db2", level=None, mode="symmetric"):
    """Single numerical backend used by both modern and compatibility interfaces."""
    if transform not in ("fourier", "z", "mellin", "wavelet"):
        raise ValueError("transform must be fourier, z, mellin, or wavelet")
    validate_preprocessing(normalization, center)
    length = _length(length, sequence_length)
    signals = {name: prepare_signal(x, normalization, center) for name, x in vectors.items()}
    metadata = dict(transform=transform, length=length, sequence_length=sequence_length,
                    normalization=normalization, center=bool(center))
    coefficients = {}
    frequency = fftshift(fftfreq(length))
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        if transform in ("fourier", "z"):
            if transform == "z":
                if isinstance(radius, (bool, complex)) or not np.isscalar(radius):
                    raise ValueError("radius must be a positive finite real scalar")
                radius = float(radius)
                if not np.isfinite(radius) or radius <= 0:
                    raise ValueError("radius must be positive and finite")
                weights = np.exp(-np.arange(sequence_length)*np.log(radius))
                metadata['radius'] = radius
            else:
                weights = 1.0
            coefficients = {name: fftshift(fft(x*weights, n=length))
                            for name, x in signals.items()}
            coordinates = {'frequency': frequency}
            if transform == 'z':
                coordinates['z'] = radius*np.exp(2j*np.pi*frequency)
        elif transform == "mellin":
            if isinstance(mellin_c, (bool, complex)) or not np.isscalar(mellin_c):
                raise ValueError("mellin_c must be a finite real scalar")
            mellin_c = float(mellin_c)
            if not np.isfinite(mellin_c):
                raise ValueError("mellin_c must be finite")
            # Grid is frozen by padded length, not each test sequence's length.
            spacing = np.log(length)/(length-1) if length > 1 else 1.0
            nu = fftshift(fftfreq(length, d=spacing))
            positions = np.arange(1, sequence_length+1, dtype=float)
            # Evaluate one frequency at a time to avoid a length-by-length kernel.
            for name, x in signals.items():
                weighted = x*positions**(mellin_c-1)
                coefficients[name] = np.asarray([
                    np.sum(weighted*np.exp(-2j*np.pi*f*np.log(positions))) for f in nu])
            coordinates = {'log_frequency': nu}
            metadata.update(mellin_c=mellin_c, position_origin=1,
                            convention='unit-spacing finite Mellin sum, negative imaginary exponent')
        else:
            try:
                import pywt
            except ImportError as error:
                raise ImportError("wavelets require PyWavelets; install Protspect[wavelets]") from error
            w = pywt.Wavelet(wavelet)
            if mode not in pywt.Modes.modes:
                raise ValueError("unsupported wavelet boundary mode")
            maximum = pywt.dwt_max_level(length, w.dec_len)
            if level is None:
                level = maximum
            if isinstance(level, bool) or not isinstance(level, Integral) or not 0 <= level <= maximum:
                raise ValueError(f"wavelet level must be an integer from 0 to {maximum}")
            for name, x in signals.items():
                bands = pywt.wavedec(np.pad(x, (0, length-len(x))), w,
                                     mode=mode, level=int(level))
                coefficients[name] = np.concatenate(bands)
            labels = [f'A{level}']+[f'D{i}' for i in range(level, 0, -1)]
            coordinates = {'band': np.concatenate([np.repeat(label, len(band))
                                                    for label, band in zip(labels, bands)]),
                           'index': np.concatenate([np.arange(len(band)) for band in bands])}
            metadata.update(wavelet=w.name, level=int(level), mode=mode,
                            pywavelets_version=pywt.__version__, padding='right zero padding')
    if any(not np.isfinite(x).all() for x in coefficients.values()):
        raise ValueError("transform produced nonfinite coefficients; rescale or change parameters")
    return TransformResult(coefficients, coordinates, metadata)


def transform_pair(first, second, embedding=b100, *, method="fourier", transform=None, length=None, **options):
    """Transform two chains on one grid and return all ordered coefficient products.

    Options follow transform_sequence. The default grid covers both chains;
    preprocessing is applied to each chain independently. Products retain the
    transform's actual coordinates (including wavelet band/position).
    """
    transform = _resolve_method(method, transform)
    first, second = _sequences([first, second])
    length = _length(length, max(len(first), len(second)))
    a = transform_sequence(first, embedding, method=transform, length=length, **options)
    b = transform_sequence(second, embedding, method=transform, length=length, **options)
    return a.cross(b)


def _resolve_method(method, alias):
    """Resolve the former transform keyword without silently accepting conflicts."""
    if alias is not None:
        if method != 'fourier' and method != alias:
            raise ValueError('method and compatibility transform arguments conflict')
        method = alias
    if method not in ('fourier', 'z', 'mellin', 'wavelet'):
        raise ValueError('method must be fourier, z, mellin, or wavelet')
    return method
