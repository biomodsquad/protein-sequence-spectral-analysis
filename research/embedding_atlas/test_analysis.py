"""Analytical checks for the unsupervised representation comparison metrics."""

import numpy as np
import pytest
from numpy.testing import assert_allclose

from research.embedding_atlas.analyze import (
    cluster_summary, deduplicate, gram_features, projector_features,
    shape_distance, unit_profile,
)


def test_vector_projectors_encode_squared_pearson_and_ignore_sign():
    a = np.array([1., 2., 4., 3.])
    b = np.array([2., 1., 3., 5.])
    u, v = unit_profile(a), unit_profile(b)
    features = projector_features(np.array([u, v, -u]))
    similarity, distances = shape_distance(features)
    r = np.corrcoef(a, b)[0, 1]
    assert_allclose(similarity[0, 1], r*r)
    assert_allclose(distances[0, 1], np.sqrt(2*(1-r*r)))
    assert_allclose(features[0], features[2])
    assert_allclose(unit_profile(a*3+8), u)
    assert unit_profile(np.ones(20)) is None


def test_cka_matches_independent_feature_space_formula():
    rng = np.random.default_rng(5)
    x, y = rng.normal(size=(20, 10)), rng.normal(size=(20, 20))
    gx, _ = gram_features(x, False)
    gy, _ = gram_features(y, False)
    x, y = x-x.mean(axis=0), y-y.mean(axis=0)
    expected = np.linalg.norm(x.T@y)**2/(np.linalg.norm(x.T@x)*np.linalg.norm(y.T@y))
    assert_allclose(gx@gy, expected)


def test_native_cka_is_invariant_to_orthogonal_channel_rotations():
    rng = np.random.default_rng(4)
    x = rng.normal(size=(20, 10))
    q, _ = np.linalg.qr(rng.normal(size=(10, 10)))
    a, _ = gram_features(x, False)
    b, _ = gram_features(3*x@q+np.arange(10), False)
    assert_allclose(a, b, atol=1e-14)


def test_equal_channel_cka_removes_individual_channel_scale_and_offsets():
    rng = np.random.default_rng(3)
    x = rng.normal(size=(20, 10))
    a, _ = gram_features(x, True)
    b, _ = gram_features(x*np.arange(1, 11)+np.arange(10), True)
    assert_allclose(a, b, atol=1e-14)
    assert gram_features(rng.normal(size=(20, 30)), True)[1]['rank'] == 19
    with pytest.raises(ValueError):
        gram_features(np.ones((20, 10)), True)


def test_complete_linkage_cut_has_the_stated_pairwise_guarantee():
    angles = np.deg2rad([0, 10, 20, 90])
    u = np.array([np.cos(angles), np.sin(angles)]).T
    _, distance = shape_distance(projector_features(u))
    _, _, labels = cluster_summary(distance, [.95], True)
    for group in set(labels[.95]):
        members = u[labels[.95] == group]
        assert np.min(abs(members@members.T)) >= .95-1e-12


def test_deduplication_keeps_antipodal_aliases_without_dropping_other_shapes():
    a, b = unit_profile([0., 1., 2.]), unit_profile([0., 2., 1.])
    unique, aliases = deduplicate(np.array([a, -a, b]))
    assert len(unique) == 2
    assert list(aliases) == [0, 0, 1]
