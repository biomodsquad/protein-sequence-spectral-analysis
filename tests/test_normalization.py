"""Mathematical normalization identities, padding policy, and paired-chain behavior."""

import numpy as np
import pytest
from numpy.testing import assert_allclose
from scipy.fft import fft, fftshift

from protspect import ComplexSet, SequenceFFT, SequenceSet

EMBEDDING = {'x': {'A': 2., 'C': 4.}, 'y': {'A': -1., 'C': 3.}}


def encoded(sequence='AAC', length=8, **options):
    result = SequenceFFT(sequence)
    result.embed_sequence(EMBEDDING)
    result.compute_spectra(length, **options)
    return result


@pytest.mark.parametrize('mode, divisor', [('none', 1), ('length', 3), ('sqrt_length', np.sqrt(3))])
def test_length_scaling_matches_dft_and_cross_scaling(mode, divisor):
    raw, scaled = encoded(), encoded(normalization=mode)
    raw.compute_cross_spectra()
    scaled.compute_cross_spectra()
    for name in EMBEDDING:
        assert_allclose(scaled.spectra[name], raw.spectra[name] / divisor)
    for name in raw.cross_spectra:
        assert_allclose(scaled.cross_spectra[name], raw.cross_spectra[name] / divisor**2)


def test_energy_normalization_parseval_and_cross_products():
    item = encoded(normalization='energy')
    item.compute_cross_spectra()
    for spectrum in item.spectra.values():
        assert_allclose(np.sum(abs(spectrum[:-1])**2) / item.spectrum_length, 1.)
    assert_allclose(item.cross_spectra['x_y'], item.spectra['x'] * item.spectra['y'].conj())


@pytest.mark.parametrize('mode', ['none', 'length', 'sqrt_length', 'energy'])
def test_centering_precedes_padding_and_normalization(mode):
    item = encoded('AC', 8, normalization=mode, center=True)
    expected = np.array([-1., 1.])
    divisor = {'none': 1., 'length': 2., 'sqrt_length': np.sqrt(2), 'energy': np.sqrt(2)}[mode]
    direct = fftshift(fft(expected / divisor, n=8))
    assert_allclose(item.spectra['x'], np.append(direct, direct[0]), atol=1e-12)
    assert_allclose(item.spectra['x'][item.frequency == 0], 0, atol=1e-12)
    assert_allclose(item.embeddings['x'], [2., 4.])  # Never mutate the raw embedding.


@pytest.mark.parametrize('center', [False, True])
def test_zero_energy_and_constant_channels_are_finite(center):
    item = SequenceFFT('AAA')
    item.embed_sequence({'zero': {'A': 0}, 'constant': {'A': 5}})
    item.compute_spectra(8, normalization='energy', center=center)
    assert_allclose(item.spectra['zero'], 0)
    if center:
        assert_allclose(item.spectra['constant'], 0)
    assert all(np.isfinite(s).all() for s in item.spectra.values())


@pytest.mark.parametrize('mode', ['none', 'length', 'sqrt_length', 'energy'])
def test_zero_padding_does_not_change_scaling_at_shared_frequencies(mode):
    short = encoded('AC', 4, normalization=mode)
    long = encoded('AC', 8, normalization=mode)
    for f in short.frequency:
        assert_allclose(short.spectra['x'][short.frequency == f],
                        long.spectra['x'][long.frequency == f], atol=1e-12)


def test_length_normalized_constant_dc_independent_of_residue_count():
    short, long = encoded('A', normalization='length'), encoded('AAAA', normalization='length')
    assert_allclose(short.spectra['x'][short.frequency == 0], 2)
    assert_allclose(long.spectra['x'][long.frequency == 0], 2)


def test_energy_normalization_handles_extreme_finite_scales():
    for score in (1e300, 1e-300):
        item = SequenceFFT('AA')
        item.embed_sequence({'x': {'A': score}})
        item.compute_spectra(4, normalization='energy')
        assert_allclose(np.sum(abs(item.spectra['x'][:-1])**2) / 4, 1.)


def test_collection_options_match_individual_and_paired_chain_transforms():
    dataset = SequenceSet(EMBEDDING, normalization='length', center=True)
    dataset.from_sequence_list(['AC'], spectrum_length=8)
    assert_allclose(dataset.monomers[0].spectra['x'],
                    encoded('AC', normalization='length', center=True).spectra['x'])
    pairs = ComplexSet(EMBEDDING, normalization='length')
    pairs.from_sequence_list(['A'], ['ACC'], spectrum_length=8)
    first, second = encoded('A', normalization='length'), encoded('ACC', normalization='length')
    assert_allclose(pairs.AB[0].spectra['A-x_B-y'], first.spectra['x'] * second.spectra['y'].conj())
    assert pairs.AB[0].spectrum_normalization == 'length'


@pytest.mark.parametrize('options', [{'normalization': 'other'}, {'center': 'false'}])
def test_bad_preprocessing_options_are_rejected(options):
    with pytest.raises(ValueError):
        encoded(**options)
    with pytest.raises(ValueError):
        SequenceSet(**options)
    with pytest.raises(ValueError):
        ComplexSet(**options)


def test_centered_constant_float_channel_is_exactly_zero_before_energy_scaling():
    item = SequenceFFT('A' * 7)
    item.embed_sequence({'x': {'A': .1}})
    item.compute_spectra(8, normalization='energy', center=True)
    assert_allclose(item.spectra['x'], 0, atol=0)
