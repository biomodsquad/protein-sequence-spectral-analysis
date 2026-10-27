"""Mathematical identities and fixed-schema checks for alternative transforms."""

import numpy as np
import pytest
from numpy.testing import assert_allclose

from protspect import SequenceFFT, transform_sequence

EMBEDDING = {'x': {'A': 1., 'R': -2., 'N': 3.}}


def test_fourier_legacy_and_z_identity():
    legacy = SequenceFFT('ARN')
    legacy.embed_sequence(EMBEDDING)
    legacy.compute_spectra(8, center=True, normalization='energy')
    for method in ('fourier', 'z'):
        result = transform_sequence('ARN', EMBEDDING, transform=method,
                                    length=8, center=True, normalization='energy')
        assert_allclose(result.coefficients['x'], legacy.spectra['x'][:-1])
        assert len(result.coordinates['frequency']) == 8


def test_z_matches_direct_definition_and_impulse():
    result = transform_sequence('ARN', EMBEDDING, transform='z', radius=1.3, length=8)
    z = result.coordinates['z']
    assert_allclose(result.coefficients['x'], 1-2/z+3/z**2, atol=1e-14)
    impulse = transform_sequence('A', EMBEDDING, transform='z', radius=.9, length=8)
    assert_allclose(impulse.coefficients['x'], 1)


def test_mellin_definition_dc_and_singleton():
    result = transform_sequence('ARN', EMBEDDING, transform='mellin', mellin_c=.5, length=8)
    nu = result.coordinates['log_frequency']
    expected = [sum(x*t**(-.5-2j*np.pi*f) for t,x in enumerate([1,-2,3],1)) for f in nu]
    assert_allclose(result.coefficients['x'], expected)
    assert_allclose(result.coefficients['x'][nu == 0], 1-2/np.sqrt(2)+3/np.sqrt(3))
    assert_allclose(transform_sequence('A', EMBEDDING, transform='mellin').coefficients['x'], [1])


def test_wavelet_reconstruction_and_coordinates():
    pywt = pytest.importorskip('pywt')
    result = transform_sequence('ARN', EMBEDDING, transform='wavelet', length=8,
                                wavelet='haar', level=3, mode='periodization')
    values = result.coefficients['x']
    bands = result.coordinates['band']
    pieces = [values[bands == band] for band in ('A3','D3','D2','D1')]
    assert_allclose(pywt.waverec(pieces, 'haar', mode='periodization'),
                    [1,-2,3,0,0,0,0,0], atol=1e-14)
    assert_allclose(np.sum(values**2), 14)
    assert_allclose(result.features('real')[0], values)
    assert len(result.features()[1]) == len(values)


@pytest.mark.parametrize('options', [dict(transform='invalid'),dict(length=2),
    dict(transform='z',radius=0),dict(transform='z',radius=np.nan),
    dict(transform='mellin',mellin_c=np.inf),dict(normalization='bad'),
    dict(transform='z',radius=1e-300),dict(center='yes')])
def test_invalid_or_nonfinite_results_rejected(options):
    with pytest.raises(ValueError):
        transform_sequence('ARN',EMBEDDING,**options)


def test_wavelet_optional_dependency(monkeypatch):
    import sys
    monkeypatch.setitem(sys.modules, 'pywt', None)
    with pytest.raises(ImportError, match='Protspect\\[wavelets\\]'):
        transform_sequence('ARN', EMBEDDING, transform='wavelet')
    assert transform_sequence('ARN', EMBEDDING).coefficients


@pytest.mark.parametrize('method', ['fourier','z','mellin','wavelet'])
def test_estimator_fixed_schema_clone_and_frozen_options(method):
    pytest.importorskip('sklearn')
    if method == 'wavelet':
        pytest.importorskip('pywt')
    from sklearn.base import clone
    from protspect.sklearn import SequenceTransformer
    options = {'radius':1.2} if method == 'z' else {}
    estimator = SequenceTransformer(embedding=EMBEDDING, method=method,
                                    length=8, transform_options=options, component='real')
    values = estimator.fit_transform(iter(['ARN','AR']))
    assert values.shape == (2,len(estimator.get_feature_names_out()))
    assert_allclose(clone(estimator).fit_transform(['ARN','AR']), values)
    options['radius'] = 2
    assert_allclose(estimator.transform(['ARN','AR']), values)
    with pytest.raises(ValueError):
        estimator.transform(['A'*9])
    estimator.set_params(length=16)
    assert not hasattr(estimator, 'feature_names_out_')


def test_estimator_channel_statistics_and_failed_refit():
    pytest.importorskip('sklearn')
    from protspect.sklearn import SequenceTransformer
    e = SequenceTransformer(embedding=EMBEDDING, standardize_channels=True).fit(['AR'])
    assert_allclose(e.channel_mean_['x'], -.5)
    before = e.channel_mean_.copy()
    e.transform(['N'])
    assert e.channel_mean_ == before
    with pytest.raises(ValueError):
        e.fit(['X'])
    assert not hasattr(e, 'settings_')


@pytest.mark.parametrize('kwargs', [dict(level=-1),dict(level=99),dict(level=True),
                                    dict(mode='invalid')])
def test_wavelet_configuration_rejected(kwargs):
    pytest.importorskip('pywt')
    with pytest.raises(ValueError):
        transform_sequence('ARN', EMBEDDING, transform='wavelet', length=8, **kwargs)


def test_wavelet_zero_level_and_zero_energy():
    pytest.importorskip('pywt')
    result = transform_sequence('AAA', EMBEDDING, transform='wavelet', length=8,
                                level=0, center=True, normalization='energy')
    assert_allclose(result.coefficients['x'], np.zeros(8))
    assert set(result.coordinates['band']) == {'A0'}
