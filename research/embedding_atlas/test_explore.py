"""Notebook controls must change the analysis and never reuse incompatible layouts."""

import sys
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from numpy.testing import assert_allclose

from research.embedding_atlas.analyze import gram_features,unit_profile
from research.embedding_atlas.explore import (
    cached_umap,cluster_profiles,prepare_vectors,prepare_whole,project_pca,umap_signature,
)


def source():
    a = unit_profile(np.arange(20,dtype=float))
    b = unit_profile(np.sin(np.arange(20)))
    c = unit_profile(np.cos(np.arange(20)))
    table = pd.DataFrame(dict(id=['a','negative-a','b','c'],family=['AAindex1']*4,
                              axis=['']*4,description=['synthetic']*4))
    return table,np.array([a,-a,b,c])


def test_sign_toggle_changes_aliases_and_metric():
    table,profiles = source()
    aliases,points,_,_,similarity,_ = prepare_vectors(table,profiles,['AAindex1'])
    assert len(points) == 3 and list(aliases.unique_id) == [0,0,1,2]
    signed = prepare_vectors(table,profiles,['AAindex1'],sign_invariant=False)
    assert len(signed[1]) == 4
    assert_allclose(signed[4][0,1],-1)
    assert_allclose(similarity[0,1],np.corrcoef(profiles[0],profiles[2])[0,1]**2)


def test_family_and_orientation_selection_rebuilds_the_candidate_pool():
    table,profiles = source()
    table.loc[3,'family'] = 'Kidera'
    with pytest.raises(ValueError,match='unique shapes'):
        prepare_vectors(table,profiles,['AAindex1'])
    table.loc[3,'axis'] = 'columns'
    with pytest.raises(ValueError):
        prepare_vectors(table,profiles,['AAindex1','Kidera'],axis='rows')
    result = prepare_vectors(table,profiles,['AAindex1','Kidera'],deduplicate=False)
    assert len(result[1]) == 4


def test_whole_weighting_control_matches_reference_cka_helper():
    rng = np.random.default_rng(8)
    reps = [dict(id=str(i),family='AAindex2',orientation='rows',square=True,
                 description='synthetic',embedding=rng.normal(size=(20,20))) for i in range(3)]
    equal = prepare_whole(reps,['AAindex2'],equal_channels=True)
    native = prepare_whole(reps,['AAindex2'],equal_channels=False)
    assert_allclose(equal[1][0],gram_features(reps[0]['embedding'],True)[0])
    assert not np.allclose(equal[2],native[2])
    with pytest.raises(ValueError,match='at least three'):
        prepare_whole(reps,['Kidera'])


def test_linkage_controls_report_only_valid_pairwise_guarantees():
    _,profiles = source()
    distance = np.linalg.norm(profiles[:,None,:]-profiles[None,:,:],axis=2)
    _,complete,_ = cluster_profiles(distance,[.8,.95],sign_invariant=False)
    _,average,_ = cluster_profiles(distance,[.8,.95],sign_invariant=False,method='average')
    assert complete.pairwise_guarantee.all() and not average.pairwise_guarantee.any()
    with pytest.raises(ValueError):
        cluster_profiles(distance,[1.1],sign_invariant=True)


def test_pca_component_control_is_bounded_by_the_selected_data():
    rng = np.random.default_rng(2)
    scores,model = project_pca(rng.normal(size=(6,20)),10)
    assert scores.shape == (6,5) and len(model.explained_variance_ratio_) == 5
    with pytest.raises(ValueError):
        project_pca(rng.normal(size=(6,20)),1)


def test_cache_key_changes_for_parameters_or_input_and_random_fits_are_not_reused(monkeypatch):
    calls = []
    class FakeUMAP:
        def __init__(self,**parameters):
            self.parameters = parameters
        def fit_transform(self,data):
            calls.append(self.parameters)
            return np.zeros((len(data),2))+len(calls)
    monkeypatch.setitem(sys.modules,'umap',SimpleNamespace(UMAP=FakeUMAP))
    distance = np.ones((5,5))-np.eye(5)
    settings = dict(n_neighbors=15,min_dist=.15,random_state=0)
    cache = {}
    first = cached_umap(distance,settings,cache)
    second = cached_umap(distance,settings,cache)
    assert not first[3] and second[3] and len(calls) == 1
    assert first[1]['n_neighbors'] == 4
    assert_allclose(first[0],second[0])
    changed = cached_umap(distance,settings | {'min_dist':.4},cache)
    assert changed[2] != first[2] and len(calls) == 2
    altered = distance.copy()
    altered[0,1] = altered[1,0] = .5
    assert umap_signature(altered,settings)[2] != first[2]
    cached_umap(altered,settings,cache)
    for _ in range(2):
        assert not cached_umap(distance,settings | {'random_state':None},cache)[3]
    assert len(calls) == 5
