"""Cross-interface mathematical parity and compatibility for the unified API."""

import json
import pickle

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal

from protspect import (ComplexSet, Sequence, SequenceFFT, SequenceSet,
                       transform_pair, transform_sequence)

EMBEDDING = {'x': {'A': 1., 'C': -2.}, 'y': {'A': 3., 'C': 0.5}}
OPTIONS = {'fourier': {}, 'z': {'radius': 1.1}, 'mellin': {'mellin_c': .5},
           'wavelet': {'wavelet': 'haar', 'level': 2, 'mode': 'periodization'}}


def check_method(method):
    if method == 'wavelet':
        pytest.importorskip('pywt')


@pytest.mark.parametrize('method', OPTIONS)
def test_single_sequence_set_and_estimator_share_coefficients_and_names(method):
    check_method(method)
    result = transform_sequence('AC', EMBEDDING, method=method, length=8, **OPTIONS[method])
    item = Sequence('AC')
    item.embed_sequence(EMBEDDING)
    item.compute_transform(8, method=method, **OPTIONS[method])
    kinds = ('coefficients', 'cross_coefficients')
    expected, names = result.features(('real', 'imaginary'), kinds=kinds)
    assert_allclose(item.features(('real', 'imaginary'), kinds=kinds)[0], expected)
    dataset = SequenceSet(EMBEDDING, method=method, transform_options=OPTIONS[method])
    dataset.from_sequence_list(iter(['AC', 'CAA']), length=8, attributes={'label': [1, 0]})
    table = dataset.export_features(('real', 'imaginary'), kinds=kinds, attributes=['label'])
    assert list(table.label) == [1, 0]
    assert_array_equal(table.columns[1:], names)
    assert_allclose(table.iloc[0, 1:].to_numpy(dtype=float), expected)
    assert dataset.length == 8
    pytest.importorskip('sklearn')
    from protspect.sklearn import SequenceTransformer
    t = SequenceTransformer(embedding=EMBEDDING, method=method, length=8,
        transform_options=OPTIONS[method], component=('real', 'imaginary'), kinds=kinds)
    actual = t.fit_transform(iter(['AC', 'CAA']))
    assert_allclose(actual, table.drop(columns='label'))
    assert_array_equal(t.get_feature_names_out(), names)
    assert json.loads(json.dumps(t.get_feature_metadata())) == t.get_feature_metadata()
    assert_allclose(pickle.loads(pickle.dumps(t)).transform(['AC', 'CAA']), actual)


@pytest.mark.parametrize('method', OPTIONS)
def test_paired_function_complex_set_and_estimator_share_ordered_products(method):
    check_method(method)
    options = OPTIONS[method]
    first = transform_sequence('AC', EMBEDDING, method=method, length=8, **options)
    second = transform_sequence('CAA', EMBEDDING, method=method, length=8, **options)
    result = transform_pair('AC', 'CAA', EMBEDDING, method=method, length=8, **options)
    assert_allclose(result.coefficients['A-x_B-y'],
                    first.coefficients['x'] * second.coefficients['y'].conj())
    reversed_pair = transform_pair('CAA', 'AC', EMBEDDING, method=method, length=8, **options)
    assert_allclose(result.coefficients['A-x_B-y'], reversed_pair.coefficients['A-y_B-x'].conj())
    dataset = ComplexSet(EMBEDDING, method=method, transform_options=options)
    dataset.from_sequence_list(['AC'], ['CAA'], length=8, attributes={'pair_id': ['p']})
    table = dataset.export_features(('real', 'imaginary'), attributes=['pair_id'])
    expected, names = result.features(('real', 'imaginary'))
    assert_allclose(table.iloc[0, 1:].to_numpy(dtype=float), expected)
    assert_array_equal(table.columns[1:], names)
    for chain, reference in [('A', first), ('B', second)]:
        assert_allclose(dataset.export_features('real', chain_ID=chain).to_numpy()[0],
                        reference.features('real')[0])
    pytest.importorskip('sklearn')
    from protspect.sklearn import SequenceTransformer
    t = SequenceTransformer(embedding=EMBEDDING, method=method, length=8, paired=True,
                            component=('real', 'imaginary'), transform_options=options)
    assert_allclose(t.fit_transform(iter([('AC', 'CAA')]))[0], expected)
    assert_array_equal(t.get_feature_names_out(), names)


def test_coordinate_selectors_preserve_original_indices_and_band_semantics():
    pytest.importorskip('pywt')
    result = transform_sequence('AC', EMBEDDING, method='wavelet', length=8,
                                wavelet='haar', level=2)
    mask = result.coordinates['band'] == 'D1'
    values, names = result.features('real', selector=lambda c: c['band'] == 'D1')
    assert_allclose(values, np.concatenate([v[mask] for v in result.coefficients.values()]))
    assert names[0].endswith(':'+str(np.flatnonzero(mask)[0]))
    original = result.coordinates['band'].copy()
    def mutate_coordinates(c):
        c['band'][:] = 'changed'
        return np.ones(len(c['band']), dtype=bool)
    result.features(selector=mutate_coordinates)
    assert_array_equal(result.coordinates['band'], original)
    for selector in ([True], lambda c: np.ones(len(c['band']))):
        with pytest.raises(ValueError, match='boolean mask'):
            result.features(selector=selector)


def test_cross_rejects_different_grids_bases_and_ambiguous_names():
    a = transform_sequence('AC', EMBEDDING, length=8)
    for b in [transform_sequence('AC', EMBEDDING, length=4),
              transform_sequence('AC', EMBEDDING, method='z', length=8),
              transform_sequence('AC', EMBEDDING, length=8, center=True)]:
        with pytest.raises(ValueError, match='matching'):
            a.cross(b)
    ambiguous = {'a': {'A': 1}, 'a_b': {'A': 1}, 'b_c': {'A': 1}, 'c': {'A': 1}}
    with pytest.raises(ValueError, match='ambiguous'):
        transform_sequence('A', ambiguous).cross()


def test_legacy_closed_fourier_arrays_and_exports_use_same_result_engine():
    item = SequenceFFT('AC')
    item.embed_sequence(EMBEDDING)
    item.compute_spectra(8)
    reference = transform_sequence('AC', EMBEDDING, length=8)
    assert_allclose(item.spectra['x'][:-1], reference.coefficients['x'])
    assert_allclose(item.spectra['x'][-1], reference.coefficients['x'][0])
    item.compute_cross_spectra()
    direct = reference.cross()
    for key in direct.coefficients:
        assert_allclose(item.cross_spectra[key][:-1], direct.coefficients[key])
    old = item.compute_spectral_properties(('spectra', 'cross_spectra'),
        select_method=lambda f: f >= 0, compute_phase=True)
    assert 'Amp_x_0.5' in old and 'Phase_x_y_0.5' in old
    assert item.result.coordinates['frequency'].size == 8
    assert item.frequency.size == 9


def test_stale_sequence_and_collection_results_invalidated():
    item = Sequence('AC')
    item.embed_sequence(EMBEDDING)
    item.compute_transform(8)
    item.embed_sequence(EMBEDDING)
    with pytest.raises(ValueError):
        item.features()
    for dataset in [SequenceSet(EMBEDDING), ComplexSet(EMBEDDING)]:
        if isinstance(dataset, ComplexSet):
            dataset.from_sequence_list(['AC'], ['A'])
            with pytest.raises(ValueError):
                dataset.from_sequence_list(['AC'], [])
        else:
            dataset.from_sequence_list(['AC'])
            with pytest.raises(ValueError):
                dataset.from_sequence_list([])
        with pytest.raises(ValueError):
            dataset.export_features()
    with pytest.raises(ValueError, match='not both'):
        SequenceSet().from_sequence_list(['A'], spectrum_length=8, length=8)


@pytest.mark.parametrize('method', ['z', 'mellin', 'wavelet'])
def test_legacy_frequency_export_is_not_applied_to_other_axes(method):
    check_method(method)
    dataset = SequenceSet(EMBEDDING, method=method, transform_options=OPTIONS[method])
    dataset.from_sequence_list(['AC'], length=8)
    assert dataset.export_features().shape[0] == 1
    with pytest.raises(ValueError, match='use export_features'):
        dataset.export_datatable()


def test_method_spelling_is_shared_and_previous_transform_keyword_still_works():
    assert_allclose(transform_sequence('AC', EMBEDDING, method='z').coefficients['x'],
                    transform_sequence('AC', EMBEDDING, transform='z').coefficients['x'])
    with pytest.raises(ValueError, match='conflict'):
        transform_sequence('AC', EMBEDDING, method='mellin', transform='z')


def test_frozen_generic_selector_options_and_paired_training_moments():
    pytest.importorskip('sklearn')
    from protspect.sklearn import SequenceTransformer
    state = [True]
    options = {'radius': 1.1}
    def selector(c):
        return c['frequency'] >= 0 if state[0] else c['frequency'] < 0
    t = SequenceTransformer(embedding=EMBEDDING, method='z', length=8, paired=True,
        selector=selector, transform_options=options, standardize_channels=True)
    samples = [('A', 'ACC'), ('AA', 'C')]
    expected = t.fit_transform(samples)
    assert_allclose(t.channel_mean_['x'], (4*1+3*(-2))/7)
    assert_allclose(t.channel_scale_['x'], np.std([1,1,-2,-2,1,1,-2]))
    state[0] = False
    options['radius'] = 2
    t.method = 'mellin'
    assert_allclose(t.transform(samples), expected)
    assert t.get_feature_metadata()['coordinates']['z']['real']
    with pytest.raises(ValueError):
        t.transform([('A', 'C', 'A')])
    with pytest.raises(ValueError):
        t.fit([('X', 'A')])
    assert not hasattr(t, 'settings_')


def test_generic_pipeline_cross_validation_fits_only_training_moments():
    pytest.importorskip('sklearn')
    from sklearn.model_selection import StratifiedKFold, cross_validate
    from sklearn.pipeline import make_pipeline
    from sklearn.linear_model import LogisticRegression
    from protspect.sklearn import SequenceTransformer, SpectralTransformer
    assert issubclass(SpectralTransformer, SequenceTransformer)
    X = ['A','AA','AAA','AAAA','C','CC','CCC','CCCC']
    y = [0,0,0,0,1,1,1,1]
    folds = list(StratifiedKFold(2, shuffle=True, random_state=7).split(X, y))
    pipeline = make_pipeline(SequenceTransformer(embedding=EMBEDDING, length=8,
        standardize_channels=True, component=('amplitude','phase'),
        kinds=('coefficients','cross_coefficients')), LogisticRegression())
    result = cross_validate(pipeline, X, y, cv=folds, return_estimator=True, error_score='raise')
    for (indices, _), model in zip(folds, result['estimator']):
        x = [EMBEDDING['x'][aa] for i in indices for aa in X[i]]
        assert_allclose(model.steps[0][1].channel_mean_['x'], np.mean(x))


@pytest.mark.parametrize('options', [dict(kinds=()), dict(component=()), dict(paired=1),
    dict(selector=[False]*8), dict(component=('real','real')), dict(kinds=('bad',))])
def test_generic_fit_rejects_invalid_schema(options):
    pytest.importorskip('sklearn')
    from protspect.sklearn import SequenceTransformer
    t = SequenceTransformer(embedding=EMBEDDING, length=8, **options)
    with pytest.raises(ValueError):
        t.fit(['AC'])
    assert not hasattr(t, 'settings_')


def test_pre_unification_fitted_spectral_pickle_is_upgraded():
    pytest.importorskip('sklearn')
    from protspect.sklearn import SpectralTransformer
    t = SpectralTransformer(embedding=EMBEDDING, spectrum_length=8,
        spectra_types=('spectra','cross_spectra'), compute_phase=True).fit(['AC','CAA'])
    expected = t.transform(['AC','CAA'])
    metadata = t.get_feature_metadata()
    # These attributes were introduced by the shared implementation.
    for name in ['settings_', 'length_', 'component_', 'kinds_', 'coordinates_',
                 'paired_', 'legacy_', 'transform_metadata_']:
        delattr(t, name)
    restored = pickle.loads(pickle.dumps(t))
    assert_allclose(restored.transform(['AC','CAA']), expected)
    assert restored.get_feature_metadata() == metadata


def test_pre_unification_generic_pickle_is_upgraded():
    pytest.importorskip('sklearn')
    from protspect.sklearn import SequenceTransformer
    t = SequenceTransformer(embedding=EMBEDDING, method='z', length=8).fit(['AC'])
    expected = t.transform(['AC'])
    for name in ['selector', 'kinds', 'paired', 'selected_mask_', 'kinds_', 'legacy_', 'paired_', 'package_version_',
                 'standardize_channels_', 'normalization_', 'center_']:
        delattr(t, name)
    restored = pickle.loads(pickle.dumps(t))
    assert_allclose(restored.transform(['AC']), expected)
    json.dumps(restored.get_feature_metadata())
    from sklearn.base import clone
    assert_allclose(clone(restored).fit_transform(['AC']), expected)


def test_component_magnitude_overflow_rejected_by_shared_exporter():
    from protspect import TransformResult
    result = TransformResult({'x': np.array([1.7e308+1.7e308j])},
        {'frequency': np.array([0.])}, {'transform': 'fourier'})
    assert np.isfinite(result.features('real')[0]).all()
    with pytest.raises(ValueError, match='nonfinite'):
        result.features('amplitude')
