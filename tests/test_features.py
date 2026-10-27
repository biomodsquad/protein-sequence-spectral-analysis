"""Scientific and API regression tests for spectral sequence features."""

import numpy as np
import pytest
from numpy.testing import assert_allclose

from protspect import ComplexSet, SequenceFFT, SequenceSet
from protspect.embeddings import b100
from protspect.utility import (
    frequency_select_physical,
    frequency_select_positive,
    frequency_selelect_positive,
)

EMBEDDING = {"score": {"A": 1.0, "C": -1.0}}


def encoded(sequence, length):
    item = SequenceFFT(sequence)
    item.embed_sequence(EMBEDDING)
    item.compute_spectra(length)
    return item


@pytest.mark.parametrize("length", [1, 2, 5, 8])
def test_fft_matches_direct_dft_and_bin_frequencies(length):
    sequence = "ACACA"[:length]
    item = encoded(sequence, length)
    signal = np.zeros(length)
    signal[:len(sequence)] = [EMBEDDING["score"][aa] for aa in sequence]
    direct = np.exp(-2j * np.pi * np.outer(item.frequency, np.arange(length))) @ signal
    assert_allclose(item.spectra["score"], direct, atol=1e-12)
    assert_allclose(item.spectra["score"][-1], item.spectra["score"][0])
    assert_allclose(item.spectra["score"][item.frequency == 0], signal.sum())
    assert_allclose(np.sum(abs(item.spectra["score"][:-1])**2) / length,
                    np.sum(signal**2))


def test_constant_and_alternating_signals_have_correct_frequency_labels():
    constant = encoded("AAAA", 4)
    alternating = encoded("ACAC", 4)
    assert constant.compute_spectral_properties()["Amp_score_0.0"] == 4
    assert constant.compute_spectral_properties()["Amp_score_0.5"] == 0
    assert alternating.compute_spectral_properties()["Amp_score_0.5"] == 4
    assert alternating.compute_spectral_properties()["Amp_score_0.0"] == 0


def test_cross_spectra_and_phase():
    item = SequenceFFT("AC")
    item.embed_sequence({"x": {"A": 1, "C": 0}, "y": {"A": 0, "C": 1}})
    item.compute_spectra(4)
    item.compute_cross_spectra()
    assert set(item.cross_spectra) == {"x_x", "x_y", "y_y"}
    assert_allclose(item.cross_spectra["x_y"],
                    item.spectra["x"] * item.spectra["y"].conj())
    props = item.compute_spectral_properties(("cross_spectra",), compute_phase=True)
    assert_allclose(props["Phase_x_y_0.25"], np.pi / 2)


def test_sequence_set_generator_metadata_and_export():
    dataset = SequenceSet(EMBEDDING)
    dataset.from_sequence_list(iter(["AC", "AAA"]), attributes={"label": [1, 2]})
    assert dataset.spectrum_length == 4
    table = dataset.export_datatable(attributes=["label"])
    assert table.shape == (2, 4)
    assert list(table.label) == [1, 2]
    assert list(table["Amp_score_0.0"]) == [0, 3]
    dataset.from_sequence_list(["A"], attributes={})
    assert len(dataset.monomers) == 1


def test_paired_cross_spectra_and_chain_exports():
    dataset = ComplexSet(EMBEDDING)
    dataset.from_sequence_list(["AC"], ["AAA"], attributes={"label": ["pair"]})
    expected = dataset.A[0].spectra["score"] * dataset.B[0].spectra["score"].conj()
    assert_allclose(dataset.AB[0].spectra["A-score_B-score"], expected)
    for chain in ("A", "B", "AB"):
        assert dataset.export_datatable(chain, attributes=["label"]).iloc[0].label == "pair"
    with pytest.raises(ValueError):
        dataset.export_datatable("other")


@pytest.mark.parametrize("sequences", [[], [""], [None], ["A", ""], "AC"])
def test_invalid_sequence_sets(sequences):
    with pytest.raises((ValueError, TypeError)):
        SequenceSet().from_sequence_list(sequences)


@pytest.mark.parametrize("length", [0, 1, -2, True, 3.5])
def test_invalid_fft_lengths(length):
    with pytest.raises(ValueError):
        encoded("AC", length)


def test_metadata_and_pairs_do_not_silently_truncate():
    with pytest.raises(ValueError):
        SequenceSet().from_sequence_list(["A", "C"], attributes={"label": [1]})
    with pytest.raises(ValueError):
        ComplexSet().from_sequence_list(["A", "C"], ["A"])
    for name in ("sequence", "spectra", "compute_spectra", "_private"):
        with pytest.raises(ValueError):
            SequenceFFT("A", {name: "invalid"})


def test_invalid_embedding_and_selector():
    item = SequenceFFT("AX")
    with pytest.raises(ValueError, match="missing residue"):
        item.embed_sequence(b100)
    with pytest.raises(ValueError, match="finite"):
        SequenceFFT("A").embed_sequence({"x": {"A": np.nan}})
    item = encoded("AC", 4)
    with pytest.raises(ValueError, match="boolean mask"):
        item.compute_spectral_properties(select_method=lambda f: [0, 1])
    with pytest.raises(ValueError, match="unavailable"):
        item.compute_spectral_properties(("cross_spectra",))
    with pytest.raises(NotImplementedError):
        ComplexSet(compute_cross_cross_spectra=True)


def test_call_order_and_stale_state():
    item = SequenceFFT("A")
    with pytest.raises(ValueError):
        item.compute_spectra(4)
    with pytest.raises(ValueError):
        item.compute_cross_spectra()
    for dataset in (SequenceSet(), ComplexSet()):
        with pytest.raises(ValueError):
            dataset.export_datatable()
    item = encoded("A", 4)
    item.compute_cross_spectra()
    item.compute_spectra(8)
    assert not hasattr(item, "cross_spectra")
    item.embed_sequence(EMBEDDING)
    assert not hasattr(item, "spectra")


def test_dyadic_selector_and_legacy_alias():
    frequencies = np.linspace(-0.5, 0.5, 9)
    assert_allclose(frequencies[frequency_select_physical(frequencies)], [0, .125, .25, .5])
    assert_allclose(frequency_select_positive(frequencies),
                    frequency_selelect_positive(frequencies))
    assert frequency_select_physical([]).size == 0


def test_blosum_symmetry_and_immutability():
    assert len(b100) == 20
    for aa, channel in b100.items():
        assert set(channel) == set(b100)
        for bb, score in channel.items():
            assert score == b100[bb][aa]
    assert b100["P"]["A"] == -2
    with pytest.raises(TypeError):
        b100["A"]["A"] = 0


def test_legacy_namespace():
    from FFTpro.features import SequenceFFT as LegacyFFT
    from FFTpro.embeddings import b100 as legacy_b100
    assert LegacyFFT is SequenceFFT
    assert legacy_b100 is b100
