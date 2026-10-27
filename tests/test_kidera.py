"""Kidera table provenance and numerical integration regression checks."""

import hashlib
import json
import pickle

import numpy as np
import pytest
from numpy.testing import assert_allclose

from protspect import ComplexSet, SequenceFFT, SequenceSet
from protspect.embeddings import b100, kidera, kidera_factors

ALPHABET = 'ARNDCQEGHILKMFPSTWYV'


def test_entire_table_matches_pinned_cran_reference():
    assert list(kidera) == [f'KF{i}' for i in range(1, 11)]
    assert all(list(values) == list(ALPHABET) for values in kidera.values())
    rows = [[kidera[f'KF{i}'][aa] for i in range(1, 11)] for aa in ALPHABET]
    # Row-major CRAN AAdata$kideraFactors, revision ee47d5d, all 200 values.
    digest = hashlib.sha256(json.dumps(rows, separators=(',', ':')).encode()).hexdigest()
    assert digest == '25aef145d7cdecc5ab652436baa67ef24940662145ef26295318b1ab5ac9fa58'
    assert_allclose(rows[0], [-1.56, -1.67, -.97, -.27, -.93, -.78, -.20, -.08, .21, -.48])
    assert_allclose(rows[1], [.22, 1.27, 1.37, 1.87, -1.70, .46, .92, -.39, .23, .93])


def test_reference_is_immutable_and_factory_returns_independent_copies():
    with pytest.raises(TypeError):
        kidera['KF1']['A'] = 0
    first, second = kidera_factors(), kidera_factors()
    first['KF1']['A'] = 100
    assert second['KF1']['A'] == kidera['KF1']['A'] == -1.56
    assert SequenceSet().embedding is b100


def test_factor_signals_match_direct_dft_and_cross_products():
    item = SequenceFFT('ARND')
    item.embed_sequence(kidera)
    item.compute_spectra(8)
    item.compute_cross_spectra()
    for channel in kidera:
        signal = np.zeros(8)
        signal[:4] = [kidera[channel][aa] for aa in 'ARND']
        direct = np.exp(-2j*np.pi*np.outer(item.frequency, np.arange(8))) @ signal
        assert_allclose(item.spectra[channel], direct, atol=1e-12)
    assert len(item.cross_spectra) == 55
    assert_allclose(item.cross_spectra['KF1_KF10'],
                    item.spectra['KF1']*item.spectra['KF10'].conj())


def test_sequence_and_paired_features_and_unknown_residues():
    data = SequenceSet(embedding=kidera_factors(), compute_cross_spectra=False)
    data.from_sequence_list(['ACDE', 'ACDF'], spectrum_length=8)
    assert data.export_datatable().shape == (2, 40)
    pairs = ComplexSet(embedding=kidera_factors())
    pairs.from_sequence_list(['ACDE'], ['ACDF'], spectrum_length=8)
    assert pairs.export_datatable().shape == (1, 400)
    with pytest.raises(ValueError):
        data.from_sequence_list(['AX'], spectrum_length=8)


def test_fitted_kidera_transformer_cloning_serialization_and_metadata():
    pytest.importorskip('sklearn')
    from sklearn.base import clone
    from protspect.sklearn import SpectralTransformer

    transformer = SpectralTransformer(embedding=kidera_factors(), spectrum_length=8)
    x = ['ACDE', 'ACDF']
    features = transformer.fit_transform(x)
    assert features.shape == (2, 40)
    assert_allclose(clone(transformer).fit_transform(x), features)
    assert_allclose(pickle.loads(pickle.dumps(transformer)).transform(x), features)
    assert transformer.embedding_ == kidera_factors()
    assert all('KF' in name for name in transformer.get_feature_names_out())
