"""Alternative transforms with explicit coefficient coordinates and conventions."""

from dataclasses import dataclass
from numbers import Integral

import numpy as np
from scipy.fft import fft, fftfreq, fftshift

from .embeddings import b100
from .features import SequenceFFT, _length
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

    def features(self, component="amplitude"):
        """Return channel-major flattened coefficients and unambiguous names.

        Select amplitude, phase (radians), real, or imaginary components. Real
        wavelets retain signed information with ``component='real'``.
        """
        functions = {"amplitude": np.abs, "phase": np.angle,
                     "real": np.real, "imaginary": np.imag}
        if component not in functions:
            raise ValueError(f"component must be one of {tuple(functions)}")
        values, names = [], []
        for channel, coefficients in self.coefficients.items():
            values.extend(functions[component](coefficients))
            names.extend(f"{channel}:{self.metadata['transform']}:{component}:{i}"
                         for i in range(len(coefficients)))
        return np.asarray(values, dtype=float), np.asarray(names, dtype=object)


def transform_sequence(sequence, embedding=b100, *, transform="fourier", length=None,
                       normalization="none", center=False, radius=1.0,
                       mellin_c=0.0, wavelet="db2", level=None, mode="symmetric"):
    """Embed a protein sequence and transform each channel on a fixed grid.

    Parameters
    ----------
    sequence : str
        Nonempty residue sequence; all residues must occur in every channel.
    embedding : mapping
        Channel-to-residue mappings, including AAindex and Kidera mappings.
    transform : {'fourier', 'z', 'mellin', 'wavelet'}
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
        Raw coefficients, coordinates, and provenance. No cross-spectrum or
        physical-frequency selector is implicitly applied to alternative axes.
    """
    if transform not in ("fourier", "z", "mellin", "wavelet"):
        raise ValueError("transform must be fourier, z, mellin, or wavelet")
    validate_preprocessing(normalization, center)
    item = SequenceFFT(sequence)
    item.embed_sequence(embedding)
    length = _length(length, len(sequence))
    signals = {name: prepare_signal(x, normalization, center)
               for name, x in item.embeddings.items()}
    metadata = dict(transform=transform, length=length, sequence_length=len(sequence),
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
                weights = np.exp(-np.arange(len(sequence))*np.log(radius))
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
            positions = np.arange(1, len(sequence)+1, dtype=float)
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
