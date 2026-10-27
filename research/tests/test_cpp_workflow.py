"""Research-only protocol and data provenance checks; not core runtime tests."""

from pathlib import Path
import json

import numpy as np
import pandas as pd
import pytest
from numpy.testing import assert_allclose

from research.compare_upstream import compare, fasta_records
from research.cpp_workflow import clean_test_overlap, jsonable, read_table, run


def test_overlap_exclusion_preserves_inputs_and_source_row_ids():
    train = pd.DataFrame({'Sequence': ['AA', 'CC'], 'Class': [1, 0], 'source_row': [0, 1]})
    test = pd.DataFrame({'Sequence': ['AA', 'AC', 'CA'], 'Class': [1, 1, 0], 'source_row': [3, 4, 5]})
    original = test.copy()
    clean, audit = clean_test_overlap(train, test)
    pd.testing.assert_frame_equal(test, original)
    assert clean.source_row.tolist() == [4, 5]
    assert audit[0]['test_source_row'] == 3
    assert 'Sequence' not in audit[0]
    test.loc[0, 'Class'] = 0
    with pytest.raises(ValueError, match='conflicting'):
        clean_test_overlap(train, test)


@pytest.mark.parametrize('content', [
    'Sequence,Class\n,Positive\nAC,Negative\n',
    'Sequence,Class\nAX,Positive\nAC,Negative\n',
    'Sequence,Class\nAA,Positive\nAA,Positive\nCC,Negative\n',
])
def test_bad_csv_does_not_silently_drop_or_impute_records(tmp_path, content):
    path = tmp_path/'bad.csv'
    path.write_text(content)
    with pytest.raises(ValueError):
        read_table(path)


def test_upstream_comparison_preserves_labels_and_duplicate_multiplicity(tmp_path):
    fasta = tmp_path/'official.fasta'
    fasta.write_text('>Positive_1\nAA\nCC\n>Negative_1\nGG\n>Negative_2\nGG\n')
    csv = tmp_path/'subset.csv'
    csv.write_text('Sequence,Class\nAACC,Positive\nGG,Negative\n')
    assert fasta_records(fasta)[0] == ('AACC', 'Positive')
    result = compare(csv, fasta, 'https://example.test/data')
    assert result['matching_labeled_records'] == 2
    assert result['official_only'] == 1
    assert result['official_only_classes'] == {'Negative': 1}
    assert result['same_labeled_multiset'] is False


def test_unavailable_library_statistics_serialize_as_json_null():
    result = jsonable({'brier': np.float64(np.nan), 'auc': np.float64(.8)})
    assert json.loads(json.dumps(result, allow_nan=False)) == {'brier': None, 'auc': .8}


def test_small_end_to_end_example_replays_without_mistic(tmp_path):
    train = tmp_path/'train.csv'
    test = tmp_path/'test.csv'
    rows = [('A'*(6+i)+'G', 'Negative') for i in range(8)]
    rows += [('R'*(6+i)+'K', 'Positive') for i in range(8)]
    pd.DataFrame(rows, columns=['Sequence', 'Class']).to_csv(train, index=False)
    pd.DataFrame([rows[0], ('AGAGA', 'Negative'), ('RKRKR', 'Positive')],
                 columns=['Sequence', 'Class']).to_csv(test, index=False)
    first = run(train, test, tmp_path/'first', with_mistic=False)
    second = run(train, test, tmp_path/'second', with_mistic=False)
    assert first == second
    a, b = [pd.read_csv(tmp_path/name/'test_predictions.csv') for name in ('first', 'second')]
    assert_allclose(a.rbf_svc_decision, b.rbf_svc_decision, atol=1e-12)
    assert a.source_row.tolist() == [1, 2]
    assert not set(pd.read_csv(tmp_path/'first/training_split.csv').query(
        'partition == "fit"').source_row) & set(pd.read_csv(tmp_path/'first/training_split.csv').query(
        'partition == "validation"').source_row)
    with pytest.raises(FileExistsError):
        run(train, test, tmp_path/'first', with_mistic=False)


def test_prepared_official_csvs_and_overlap_match_source_manifest():
    root = Path(__file__).resolve().parents[2]
    train_path = root/'research/data/mlcpp2/training.csv'
    test_path = root/'research/data/mlcpp2/independent.csv'
    if not train_path.exists():
        pytest.skip('prepare the official website data to run this check')
    train, test = read_table(train_path), read_table(test_path)
    clean, exclusions = clean_test_overlap(train, test)
    assert len(train) == 1146
    assert len(test) == 2341
    assert len(clean) == 2337
    assert clean.Class.sum() == 153
    assert len(exclusions) == 4


def test_official_preparation_verifies_hashes_and_preserves_identifiers(tmp_path, monkeypatch):
    import hashlib
    from research import prepare_cpp_data as preparation
    cache = tmp_path/'cache'
    cache.mkdir()
    spec = {'source_page': 'https://example.test/', 'snapshot_date': '2026-10-27', 'datasets': {}}
    for name in ('training', 'independent'):
        raw = b'>Positive_1\nAA\n>Negative_1\nCC\n'
        (cache/f'mlcpp2_{name}.fasta').write_bytes(raw)
        spec['datasets'][name] = {'url': 'https://example.test/'+name,
            'sha256': hashlib.sha256(raw).hexdigest(), 'rows': 2,
            'classes': {'Positive': 1, 'Negative': 1}}
    manifest = tmp_path/'manifest.json'
    manifest.write_text(json.dumps(spec))
    monkeypatch.setattr(preparation, 'MANIFEST', manifest)
    result = preparation.prepare(tmp_path/'prepared', cache)
    table = read_table(tmp_path/'prepared/training.csv')
    assert table.source_id.tolist() == ['Positive_1', 'Negative_1']
    assert table.Sequence.tolist() == ['AA', 'CC']
    assert result['datasets']['training']['rows'] == 2
    with pytest.raises(FileExistsError):
        preparation.prepare(tmp_path/'prepared', cache)
    (cache/'mlcpp2_training.fasta').write_text('tampered')
    with pytest.raises(ValueError, match='SHA-256'):
        preparation.prepare(tmp_path/'bad', cache)
    assert not (tmp_path/'bad').exists()


def test_official_parser_rejects_duplicate_identifiers_or_sequences():
    from research.prepare_cpp_data import parse_fasta
    for text in ('>Positive_1\nAA\n>Positive_1\nAC\n',
                 '>Positive_1\nAA\n>Negative_1\nAA\n',
                 '>Unknown_1\nAA\n'):
        with pytest.raises(ValueError):
            parse_fasta(text)
