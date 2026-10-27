"""Original synthetic matrices test orientation, reconstruction, and spectral use."""

import io
import json
import pickle

import numpy as np
import pytest
from numpy.testing import assert_allclose

from protspect import ComplexSet, SequenceFFT, SequenceSet
from protspect.aaindex import (
    AAIndex, AAIndex2, AAIndex3, AAINDEX2_URL, AAINDEX3_URL,
    download_aaindex2, download_aaindex3,
)

ALPHABET = 'ARNDCQEGHILKMFPSTWYV'


def record(*, rows='AR', columns='AR', data=('1', '2 3'), accession='TEST000001'):
    """Construct original test data, not redistributed AAindex records."""
    return (f'H {accession}\nD Synthetic matrix\nR PMID:1\nA Example, A.\n'
            'T Matrix title\n  continued\nJ Example Journal\n'
            '* See https://example.org/matrix\n'
            f'M rows = {rows}, cols = {columns}\n'
            + '\n'.join('  '+line for line in data)+'\n//\n')


def load(tmp_path, text=None, cls=AAIndex2):
    p = tmp_path/'matrix'
    p.write_text(record() if text is None else text, encoding='utf-8')
    return cls.from_file(p)


def full_record(*, extra=False, missing=False, accession='TEST000001'):
    rows = ALPHABET + ('J-' if extra else '')
    columns = ALPHABET + ('J' if extra else '')
    matrix = [[100*i+j for j in range(len(columns))] for i in range(len(rows))]
    if missing:
        matrix[0][1] = 'NA'
    return record(rows=rows, columns=columns, accession=accession,
                  data=tuple(' '.join(map(str, row)) for row in matrix))


@pytest.mark.parametrize('cls', [AAIndex2, AAIndex3])
def test_lower_triangle_reconstruction_and_citation(tmp_path, cls):
    db = load(tmp_path, cls=cls)
    e = db['TEST000001']
    assert e.shape == (2, 2) and e.storage == 'lower_triangle'
    assert e.values['A']['R'] == e.values['R']['A'] == 2
    assert e.values['A']['A'] == 1 and e.values['R']['R'] == 3
    assert e.complete and not e.canonical_complete
    assert e.title == 'Matrix title continued'
    assert db.search('SYNTHETIC') == (e,)
    assert db.metadata('TEST000001')['database'] == cls.database_name
    assert db.metadata('TEST000001')['representation'] is None
    with pytest.raises(TypeError):
        e.values['A']['R'] = 0


def test_asymmetric_rectangular_full_matrix_is_preserved(tmp_path):
    e = load(tmp_path, record(rows='RA-', columns='AR', data=('1 2', '3 4', '5 6')))['TEST000001']
    assert e.shape == (3, 2) and e.storage == 'full'
    assert e.values['R']['A'] == 1
    assert e.values['A']['R'] == 4  # Not reflected from the other direction.
    assert e.values['-']['R'] == 6


def test_canonical_profiles_preserve_direction_and_exclude_extra_symbols(tmp_path):
    db = load(tmp_path, full_record(extra=True))
    rows = db.embedding('TEST000001', axis='rows')
    cols = db.embedding('TEST000001', axis='columns')
    assert len(rows) == len(cols) == 20
    assert rows['TEST000001_rows_A']['R'] == 1
    assert cols['TEST000001_columns_A']['R'] == 100
    assert list(rows) == [f'TEST000001_rows_{aa}' for aa in ALPHABET]
    assert list(rows['TEST000001_rows_A']) == list(ALPHABET)
    rows['TEST000001_rows_A']['R'] = -1
    assert db['TEST000001'].values['A']['R'] == 1
    meta = db.metadata(iter(['TEST000001']), axis='columns')
    assert meta['channels'] == list(cols)
    assert meta['entries'][0]['rows'].endswith('J-')
    assert meta['representation']['axis'] == 'columns'
    assert meta['representation']['extra_symbols'] == 'excluded'
    json.dumps(meta, allow_nan=False)
    with pytest.raises(TypeError):
        db.embedding('TEST000001')
    with pytest.raises(ValueError, match='axis'):
        db.embedding('TEST000001', axis='invalid')


def test_profiles_follow_labels_even_when_native_axes_are_reordered(tmp_path):
    # Row A is the final native row; column A is the first native column.
    native = ALPHABET[::-1]
    text = record(rows=native, columns=ALPHABET,
                  data=tuple(' '.join(str(i*100+j) for j in range(20)) for i in range(20)))
    db = load(tmp_path, text)
    assert db.embedding('TEST000001', axis='rows')['TEST000001_rows_A']['R'] == 1901
    assert db.embedding('TEST000001', axis='columns')['TEST000001_columns_A']['R'] == 1800


@pytest.mark.parametrize('marker', ['NA', '-'])
def test_missing_markers_reconstruct_without_imputation(tmp_path, marker):
    db = load(tmp_path, record(data=('1', f'{marker} 3')), cls=AAIndex3)
    e = db['TEST000001']
    assert not e.complete
    assert e.missing_pairs == (('A', 'R'), ('R', 'A'))
    assert e.values['A']['R'] is None
    assert db.search('synthetic', complete_only=True) == ()
    assert db.metadata('TEST000001')['entries'][0]['values']['A']['R'] is None


def test_incomplete_canonical_profiles_are_rejected(tmp_path):
    db = load(tmp_path, full_record(missing=True))
    assert not db['TEST000001'].canonical_complete
    for axis in ('rows', 'columns'):
        with pytest.raises(ValueError, match='missing canonical matrix values'):
            db.embedding('TEST000001', axis=axis)
        with pytest.raises(ValueError, match='missing canonical matrix values'):
            db.metadata('TEST000001', axis=axis)
    with pytest.raises(ValueError, match='canonical residues absent'):
        load(tmp_path).embedding('TEST000001', axis='rows')


def test_missing_extra_cell_does_not_invalidate_complete_canonical_block(tmp_path):
    text = full_record(extra=True)
    # The final row contains only noncanonical '-' scores.
    text = text.replace('2100', 'NA')
    db = load(tmp_path, text)
    assert not db['TEST000001'].complete and db['TEST000001'].canonical_complete
    assert len(db.embedding('TEST000001', axis='rows')) == 20


def test_multiple_matrices_keep_selected_order_and_accession_namespaces(tmp_path):
    db = load(tmp_path, full_record()+full_record(accession='TEST000002'))
    embedding = db.embedding(['TEST000002', 'TEST000001'], axis='rows')
    assert len(embedding) == 40
    assert next(iter(embedding)) == 'TEST000002_rows_A'
    with pytest.raises(ValueError, match='duplicate'):
        db.embedding(['TEST000001', 'TEST000001'], axis='rows')
    with pytest.raises(KeyError):
        db.embedding('NONE000001', axis='rows')


@pytest.mark.parametrize('text,match', [
    (record(rows='AA'), 'duplicate symbols'),
    (record(data=('1 2',)), 'do not match'),
    (record(rows='AR', columns='RA'), 'do not match'),
    (record(data=('1 2', '3')), 'do not match'),  # Upper triangle is not lower input.
    (record(data=('1', 'inf 3')), 'nonfinite'),
    (record(data=('1', 'bad 3')), 'invalid matrix value'),
    (record().replace('rows = AR', 'rows = A1'), 'axis header'),
    (record().replace('M rows', 'I rows'), 'requires M'),
])
def test_malformed_matrix_data_is_rejected(tmp_path, text, match):
    with pytest.raises(ValueError, match=match):
        load(tmp_path, text)


def test_underscore_accession_and_scalar_matrix_separation(tmp_path):
    assert list(load(tmp_path, record(accession='QU_C930101'))) == ['QU_C930101']
    with pytest.raises(ValueError, match='matrix records'):
        load(tmp_path, cls=AAIndex)


@pytest.mark.parametrize('download,cls,url', [
    (download_aaindex2, AAIndex2, AAINDEX2_URL),
    (download_aaindex3, AAIndex3, AAINDEX3_URL),
])
def test_download_helpers_use_correct_section_and_preserve_cache(tmp_path, monkeypatch, download,
                                                                 cls, url):
    requests = []
    def retrieve(address, timeout):
        requests.append(address)
        return io.BytesIO(record().encode())
    monkeypatch.setattr('protspect.aaindex.urlopen', retrieve)
    path = tmp_path/'cache'
    db = download(path)
    assert isinstance(db, cls) and db.source_url == url and requests == [url]
    assert db.sha256 == cls.from_file(path).sha256
    with pytest.raises(FileExistsError):
        download(path)
    with pytest.raises(ValueError, match='SHA256 mismatch'):
        download(tmp_path/'bad', expected_sha256='0'*64)
    assert not (tmp_path/'bad').exists()


def test_matrix_profiles_match_direct_dft_and_feature_classes(tmp_path):
    db = load(tmp_path, full_record())
    embedding = db.embedding('TEST000001', axis='columns')
    item = SequenceFFT('AR')
    item.embed_sequence(embedding)
    item.compute_spectra(4)
    signal = np.array([0, 100, 0, 0])
    expected = np.exp(-2j*np.pi*np.outer(item.frequency, np.arange(4))) @ signal
    assert_allclose(item.spectra['TEST000001_columns_A'], expected, atol=1e-12)
    single = SequenceSet(embedding=embedding, compute_cross_spectra=False)
    single.from_sequence_list(['AR'], spectrum_length=8)
    assert single.export_datatable().shape == (1, 80)
    paired = ComplexSet(embedding=embedding)
    paired.from_sequence_list(['AR'], ['RA'], spectrum_length=8)
    assert paired.export_datatable().shape == (1, 1600)


def test_matrix_profiles_support_fitted_estimator_clone_and_pickle(tmp_path):
    pytest.importorskip('sklearn')
    from sklearn.base import clone
    from protspect.sklearn import SpectralTransformer
    db = load(tmp_path, full_record(), cls=AAIndex3)
    model = SpectralTransformer(embedding=db.embedding('TEST000001', axis='rows'),
                                spectrum_length=8, standardize_channels=True)
    x = ['AR', 'RA']
    expected = model.fit_transform(x)
    assert_allclose(clone(model).fit_transform(x), expected)
    assert_allclose(pickle.loads(pickle.dumps(model)).transform(x), expected)
