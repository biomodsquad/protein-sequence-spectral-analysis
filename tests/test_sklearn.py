"""Fitted feature schemas, preprocessing leakage, and estimator interoperability."""

import pickle

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal

pytest.importorskip('sklearn')
from sklearn.base import clone
from sklearn.exceptions import NotFittedError
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_validate, GridSearchCV
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from protspect import SequenceSet
from protspect.sklearn import SpectralTransformer

EMBEDDING = {'x': {'A': 0., 'C': 6.}}


def test_default_transformer_matches_existing_exports_and_fixed_order():
    sequences = ['ACDE', 'ACDF']
    transformer = SpectralTransformer(spectrum_length=8)
    actual = transformer.fit_transform(sequences)
    dataset = SequenceSet()
    dataset.from_sequence_list(sequences, spectrum_length=8)
    expected = dataset.export_datatable()
    assert_allclose(actual, expected.to_numpy())
    assert_array_equal(transformer.get_feature_names_out(), expected.columns)


def test_grid_and_columns_do_not_depend_on_later_batch_lengths():
    transformer = SpectralTransformer().fit(['ACD'])
    assert transformer.spectrum_length_ == 4
    names = transformer.get_feature_names_out()
    single = transformer.transform(['A'])
    batch = transformer.transform(['A', 'ACDE'])
    assert_allclose(single[0], batch[0])
    assert_array_equal(names, transformer.get_feature_names_out())
    with pytest.raises(ValueError, match='exceeds fitted'):
        transformer.transform(['AAAAA'])
    assert_array_equal(names, transformer.get_feature_names_out())


def test_channel_statistics_are_residue_weighted_and_never_learn_from_test():
    transformer = SpectralTransformer(embedding=EMBEDDING, spectrum_length=8,
                                     standardize_channels=True).fit(['A', 'AAA', 'CC'])
    assert_allclose(transformer.channel_mean_['x'], 2.)
    assert_allclose(transformer.channel_scale_['x'], np.sqrt(8))
    actual = transformer.transform(['CCCC'])
    dc = list(transformer.get_feature_names_out()).index('Amp_x_0.0')
    assert_allclose(actual[0, dc], 4 * (6. - 2.) / np.sqrt(8))
    assert_allclose(transformer.channel_mean_['x'], 2.)
    assert_allclose(transformer.channel_scale_['x'], np.sqrt(8))


def test_zero_variance_uses_unit_scale_and_handles_unseen_valid_residue():
    transformer = SpectralTransformer(embedding=EMBEDDING, spectrum_length=4,
                                     standardize_channels=True).fit(['AA'])
    assert transformer.channel_scale_['x'] == 1
    assert_allclose(transformer.transform(['AA']), 0)
    dc = list(transformer.get_feature_names_out()).index('Amp_x_0.0')
    assert_allclose(transformer.transform(['C'])[0, dc], 6)


def test_freezes_custom_embedding_selector_and_constructor_parameters():
    embedding = {'x': {'A': 1., 'C': 2.}}
    state = [True]

    def selector(frequencies):
        return frequencies >= 0 if state[0] else frequencies < 0

    transformer = SpectralTransformer(embedding=embedding, select_method=selector).fit(['AC'])
    expected = transformer.transform(['AC'])
    embedding['x']['A'] = 100.
    state[0] = False
    transformer.normalization = 'length'
    assert_allclose(transformer.transform(['AC']), expected)
    assert transformer.embedding_['x']['A'] == 1.
    names = transformer.get_feature_names_out()
    names[0] = 'edited'
    assert transformer.get_feature_names_out()[0] != 'edited'


def test_not_fitted_parameter_changes_and_failed_refits():
    transformer = SpectralTransformer()
    with pytest.raises(NotFittedError):
        transformer.transform(['AC'])
    transformer.fit(['AC']).set_params(normalization='length')
    with pytest.raises(NotFittedError):
        transformer.get_feature_names_out()
    transformer.fit(['AC'])
    with pytest.raises(ValueError):
        transformer.fit(['AZ'])
    with pytest.raises(NotFittedError):
        transformer.transform(['AC'])


def test_cloning_persistence_and_single_pass_generator_fit_transform():
    transformer = SpectralTransformer(normalization='energy', center=True)
    actual = transformer.fit_transform(iter(['ACDE', 'ACDF']))
    restored = pickle.loads(pickle.dumps(transformer))
    assert_allclose(restored.transform(['ACDE', 'ACDF']), actual)
    copied = clone(transformer)
    assert copied.get_params() == transformer.get_params()
    assert not hasattr(copied, 'spectrum_length_')


def test_pipeline_cv_fits_channel_statistics_on_each_training_fold():
    X = ['A', 'AA', 'AAA', 'AAAA', 'C', 'CC', 'CCC', 'CCCC']
    y = np.array([0, 0, 0, 0, 1, 1, 1, 1])
    folds = list(StratifiedKFold(2, shuffle=True, random_state=7).split(X, y))
    pipeline = make_pipeline(SpectralTransformer(embedding=EMBEDDING, spectrum_length=8,
        standardize_channels=True), StandardScaler(), LogisticRegression())
    result = cross_validate(pipeline, X, y, cv=folds, return_estimator=True, error_score='raise')
    for (train, _), estimator in zip(folds, result['estimator']):
        values = [EMBEDDING['x'][aa] for index in train for aa in X[index]]
        features = estimator.steps[0][1]
        assert_allclose(features.channel_mean_['x'], np.mean(values))
        assert_allclose(features.channel_scale_['x'], np.std(values))
    search = GridSearchCV(pipeline, {'spectraltransformer__normalization': ['none', 'length']},
                          cv=folds, error_score='raise').fit(X, y)
    assert search.predict(['AAA', 'CCC']).shape == (2,)


def test_pandas_output_has_fitted_feature_columns():
    transformer = SpectralTransformer(spectrum_length=4).set_output(transform='pandas')
    table = transformer.fit_transform(['AC'])
    assert_array_equal(table.columns, transformer.get_feature_names_out())


@pytest.mark.parametrize('options', [
    {'normalization': 'invalid'}, {'center': 'true'}, {'standardize_channels': 1},
    {'compute_amplitude': False}, {'spectra_types': ()}, {'spectra_types': ('other',)},
    {'spectra_types': ('spectra', 'spectra')}, {'select_method': lambda f: np.ones(len(f))},
    {'select_method': lambda f: np.zeros(len(f), dtype=bool)},
    {'embedding': {'x': {'A': np.nan}}}, {'embedding': {'x': {'A': 1j}}},
    {'embedding': {'x': {'A': [1]}}}, {'embedding': {'x': {'AA': 1}}},
    {'embedding': {'x': {'A': 1}, 'y': {'C': 1}}},
])
def test_invalid_configuration_rejected_during_fit(options):
    with pytest.raises(ValueError):
        SpectralTransformer(**options).fit(['A'])


def test_invalid_transform_sequences_are_rejected_without_changing_fit():
    transformer = SpectralTransformer().fit(['AC'])
    for X in ([], [''], ['X'], 'AC'):
        with pytest.raises((ValueError, TypeError)):
            transformer.transform(X)
    assert transformer.transform(['AC']).shape[1] == transformer.n_features_out_


def test_constant_training_channel_and_metadata_round_trip():
    import json
    transformer = SpectralTransformer(embedding={'x': {'A': .1, 'C': .1}},
        standardize_channels=True, normalization='energy', center=True,
        spectrum_length=4).fit(['A', 'CC'])
    assert transformer.channel_scale_['x'] == 1.
    assert_allclose(transformer.transform(['ACC']), 0, atol=0)
    metadata = transformer.get_feature_metadata()
    assert json.loads(json.dumps(metadata)) == metadata
    metadata['embedding']['x']['A'] = 9
    metadata['feature_names'][0] = 'changed'
    assert transformer.embedding_['x']['A'] == .1
    assert transformer.get_feature_names_out()[0] != 'changed'


@pytest.mark.parametrize('score', [1e300, 1e-300])
def test_channel_standardization_avoids_squaring_extreme_scores(score):
    transformer = SpectralTransformer(embedding={'x': {'A': score, 'C': -score}},
        standardize_channels=True).fit(['AC'])
    assert_allclose(transformer.channel_scale_['x'] / score, 1)
    assert_allclose(transformer.effective_embedding_['x']['A'], 1)
    assert_allclose(transformer.effective_embedding_['x']['C'], -1)
    assert np.isfinite(transformer.transform(['AC'])).all()
