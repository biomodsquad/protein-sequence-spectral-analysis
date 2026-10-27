"""Fourier features for individual sequences, sequence sets, and paired chains."""

from collections.abc import Mapping
from numbers import Integral

import numpy as np
import pandas as pd
from scipy.fft import fft, fftfreq, fftshift

from .embeddings import b100
from .utility import frequency_select_physical
from .preprocessing import prepare_signal, validate_preprocessing


def _sequences(sequences):
    """Materialize and validate a nonempty iterable of nonempty strings."""
    if isinstance(sequences, str):
        raise TypeError("sequences must be an iterable of strings, not a string")
    result = list(sequences)
    if not result or any(not isinstance(s, str) or not s for s in result):
        raise ValueError("provide at least one nonempty sequence string")
    return result


def _length(length, minimum):
    """Choose a power-of-two length or validate an explicit FFT length."""
    if length is None:
        return 1 << (minimum - 1).bit_length()
    if isinstance(length, bool) or not isinstance(length, Integral) or length < minimum:
        raise ValueError("spectrum length must be an integer at least as long as every sequence")
    return int(length)


def _attributes(attributes, count):
    """Convert column-oriented metadata to rows without silently truncating."""
    if attributes is None:
        return [None] * count
    if not isinstance(attributes, Mapping):
        raise TypeError("attributes must be a mapping of names to columns")
    columns = {name: list(values) for name, values in attributes.items()}
    if any(len(values) != count for values in columns.values()):
        raise ValueError("every attribute column must match the number of sequences")
    return [{name: values[i] for name, values in columns.items()} for i in range(count)]


class SequenceFFT:
    """Encode a sequence using numeric residue-property channels and an FFT.

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
                    "spectra", "cross_spectra", "spectrum_normalization", "spectrum_centered"}
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
        if not isinstance(embeddings, Mapping) or not embeddings:
            raise ValueError("embeddings must be a nonempty mapping of channels")
        vectors = {}
        for name, values in embeddings.items():
            if not isinstance(name, str) or not name:
                raise ValueError("embedding channel names must be nonempty strings")
            try:
                vector = np.asarray([values[aa] for aa in self.sequence], dtype=float)
            except KeyError as error:
                raise ValueError(f"missing residue {error.args[0]!r} in channel {name!r}") from error
            if vector.ndim != 1 or not np.all(np.isfinite(vector)):
                raise ValueError("embedding values must be finite real scalars")
            vectors[name] = vector
        self.embeddings = vectors
        for name in ("spectra", "cross_spectra", "frequency", "spectrum_length",
                     "spectrum_normalization", "spectrum_centered"):
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
        validate_preprocessing(normalization, center)
        if not hasattr(self, "embeddings"):
            raise ValueError("call embed_sequence before compute_spectra")
        length = _length(length, len(self.sequence))
        bins = fftshift(fftfreq(length))
        spectra = {}
        for name, vector in self.embeddings.items():
            signal = prepare_signal(vector, normalization, center)
            spectrum = fftshift(fft(signal, n=length))
            if not np.all(np.isfinite(spectrum)):
                raise ValueError("FFT produced nonfinite values; rescale the embedding")
            spectra[name] = np.append(spectrum, spectrum[0])
        self.spectrum_length = length
        self.spectrum_normalization = normalization
        self.spectrum_centered = bool(center)
        self.frequency = np.append(bins, bins[0] + 1)
        self.spectra = spectra
        self.__dict__.pop("cross_spectra", None)

    def compute_cross_spectra(self):
        """Store upper-triangle channel products ``X_i * conjugate(X_j)``.

        Includes auto-spectra; keys are ``channel_i_channel_j``. Products
        inherit both input channels' preprocessing. There is no additional
        normalization or statistical averaging of cross-products.
        """
        if not hasattr(self, "spectra"):
            raise ValueError("call compute_spectra before compute_cross_spectra")
        names = list(self.spectra)
        result = {}
        for i, first in enumerate(names):
            for second in names[i:]:
                key = f"{first}_{second}"
                if key in result:
                    raise ValueError(f"ambiguous cross-spectrum name: {key}")
                result[key] = self.spectra[first] * np.conjugate(self.spectra[second])
        self.cross_spectra = result

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
        selected = np.asarray(select_method(self.frequency))
        if selected.dtype != bool or selected.shape != self.frequency.shape:
            raise ValueError("select_method must return a boolean mask matching frequency")
        result = {}
        for kind in spectra_types:
            if kind not in ("spectra", "cross_spectra") or not hasattr(self, kind):
                raise ValueError(f"unavailable spectrum type: {kind!r}")
            for name, spectrum in getattr(self, kind).items():
                for enabled, prefix, values in (
                    (compute_amplitude, "Amp", np.abs(spectrum)),
                    (compute_phase, "Phase", np.angle(spectrum)),
                ):
                    if enabled:
                        for frequency, value in zip(self.frequency[selected], values[selected]):
                            key = f"{prefix}_{name}_{frequency}"
                            if key in result:
                                raise ValueError(f"duplicate feature name: {key}")
                            result[key] = value
        return result


class SequenceSet:
    """Compute features for sequences sharing one embedding and FFT grid.

    Parameters
    ----------
    embedding : mapping, optional
        Channel-to-residue score mappings; defaults to canonical BLOSUM100.
    compute_cross_spectra : bool, optional
        Compute within-sequence channel cross-products, including auto-spectra.
    normalization : {"none", "length", "sqrt_length", "energy"}, optional
        Per-channel scaling passed to SequenceFFT.compute_spectra.
    center : bool, optional
        Center each channel over actual residues before scaling and padding.
    """

    def __init__(self, embedding=b100, compute_cross_spectra=True, *,
                 normalization="none", center=False):
        self.embedding = embedding
        self.compute_cross_spectra = compute_cross_spectra
        validate_preprocessing(normalization, center)
        self.normalization = normalization
        self.center = center

    def from_sequence_list(self, sequences, spectrum_length=None, attributes=None):
        """Populate the set from an iterable of nonempty sequence strings.

        ``spectrum_length`` defaults to the smallest power of two covering the
        longest sequence. ``attributes`` maps metadata names to equal-length
        columns, one value per sequence. Mutates this object and returns None.
        """
        sequences = _sequences(sequences)
        length = _length(spectrum_length, max(map(len, sequences)))
        rows = _attributes(attributes, len(sequences))
        monomers = []
        for sequence, row in zip(sequences, rows):
            item = SequenceFFT(sequence, row)
            item.embed_sequence(self.embedding)
            item.compute_spectra(length, normalization=self.normalization, center=self.center)
            if self.compute_cross_spectra:
                item.compute_cross_spectra()
            monomers.append(item)
        self.sequences, self.spectrum_length, self.monomers = sequences, length, monomers

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
    chain independently before cross-products.
    """

    def __init__(self, embedding=b100, compute_cross_cross_spectra=False, *,
                 normalization="none", center=False):
        if compute_cross_cross_spectra:
            raise NotImplementedError("cross-cross spectra are not implemented")
        self.embedding = embedding
        self.compute_cross_cross_spectra = compute_cross_cross_spectra
        validate_preprocessing(normalization, center)
        self.normalization = normalization
        self.center = center

    def from_sequence_list(self, A_sequences, B_sequences, spectrum_length=None, attributes=None):
        """Populate equally sized collections of paired A/B sequence strings.

        All chains share the supplied FFT length, or the smallest power of two
        covering the longest chain. Attribute columns have one value per pair.
        ``A`` and ``B`` store chain spectra; ``AB`` stores every channel product
        ``FFT(A_i) * conjugate(FFT(B_j))`` after each chain's preprocessing.
        """
        a_sequences, b_sequences = _sequences(A_sequences), _sequences(B_sequences)
        if len(a_sequences) != len(b_sequences):
            raise ValueError("A_sequences and B_sequences must have equal lengths")
        length = _length(spectrum_length, max(map(len, a_sequences + b_sequences)))
        rows = _attributes(attributes, len(a_sequences))
        a_items, b_items, pairs = [], [], []
        for a, b, row in zip(a_sequences, b_sequences, rows):
            first, second = SequenceFFT(a, row), SequenceFFT(b, row)
            for item in (first, second):
                item.embed_sequence(self.embedding)
                item.compute_spectra(length, normalization=self.normalization, center=self.center)
            pair = SequenceFFT(a, row)
            pair.sequence = (a, b)
            pair.frequency, pair.spectrum_length = first.frequency, length
            pair.spectrum_normalization = self.normalization
            pair.spectrum_centered = bool(self.center)
            pair.spectra = {}
            for a_name, a_spectrum in first.spectra.items():
                for b_name, b_spectrum in second.spectra.items():
                    key = f"A-{a_name}_B-{b_name}"
                    if key in pair.spectra:
                        raise ValueError(f"ambiguous paired spectrum name: {key}")
                    pair.spectra[key] = a_spectrum * np.conjugate(b_spectrum)
            a_items.append(first)
            b_items.append(second)
            pairs.append(pair)
        self.A_sequences, self.B_sequences = a_sequences, b_sequences
        self.spectrum_length = length
        self.A, self.B, self.AB = a_items, b_items, pairs

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
        return _datatable(getattr(self, chain_ID), ("spectra",), attributes,
                          select_method, compute_amplitude, compute_phase)
