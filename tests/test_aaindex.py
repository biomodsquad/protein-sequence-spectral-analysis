"""Synthetic AAindex1 format, provenance, download, and numerical integration tests."""

import hashlib
import io
import json
import pickle

import numpy as np
import pytest
from numpy.testing import assert_allclose

from protspect import ComplexSet, SequenceFFT, SequenceSet
from protspect.aaindex import AAIndex, AAINDEX1_URL, download_aaindex1

ALPHABET = 'ARNDCQEGHILKMFPSTWYV'
HEADER = 'I A/L R/K N/M D/F C/P Q/S E/T G/W H/Y I/V'


def record(accession='TEST000001', *, values=None, header=HEADER):
    """Small original fixture; no database records are redistributed in the package."""
    values = list(range(20)) if values is None else values
    return (f'H {accession}\nD Synthetic property\nR\nA Example, A.\n'
            'T Synthetic title\n  continued title\nJ Example Journal\n'
            '* See https://example.org/property\n* Second comment\n'
            'C TEST000002 0.9\n  TEST000003 -0.8\n' + header + '\n  '
            + ' '.join(map(str, values[:10])) + '\n  '
            + ' '.join(map(str, values[10:])) + '\n//\n')


def database(tmp_path, text=None, **kwargs):
    path = tmp_path/'aaindex1'
    path.write_bytes((record() if text is None else text).encode('utf-8'))
    return AAIndex.from_file(path, **kwargs)


def test_metadata_continuations_empty_reference_comments_and_checksum(tmp_path):
    text = record()
    digest = hashlib.sha256(text.encode()).hexdigest()
    db = database(tmp_path, text, expected_sha256=digest.upper(), version='fixture',
                  source_url='https://example.org/aaindex1')
    assert len(db) == 1 and list(db) == ['TEST000001']
    entry = db['TEST000001']
    assert entry.title == 'Synthetic title continued title'
    assert entry.references == ''
    assert 'https://example.org/property' in entry.comments
    assert entry.correlations == 'TEST000002 0.9 TEST000003 -0.8'
    assert list(entry.values) == list(ALPHABET)
    assert entry.values['A'] == 0 and entry.values['L'] == 10
    assert entry.complete
    meta = json.loads(json.dumps(db.metadata('TEST000001'), allow_nan=False))
    assert meta['sha256'] == digest
    assert meta['version'] == 'fixture'
    assert meta['entries'][0]['values']['V'] == 19
    assert db.search('CONTINUED') == (entry,)
    assert db.search('example, a.') == (entry,)
    assert not db.search('absent')
    with pytest.raises(TypeError):
        entry.values['A'] = 1


def test_header_residue_order_is_decoded_instead_of_assumed(tmp_path):
    header = 'I L/A K/R M/N F/D P/C S/Q T/E W/G Y/H V/I'
    entry = database(tmp_path, record(header=header))['TEST000001']
    assert entry.values['L'] == 0 and entry.values['A'] == 10
    assert entry.values['V'] == 9 and entry.values['I'] == 19


def test_selection_preserves_order_and_returns_fresh_dicts(tmp_path):
    db = database(tmp_path, record() + record('TEST000002'))
    values = db.embedding(['TEST000002', 'TEST000001'])
    assert list(values) == ['TEST000002', 'TEST000001']
    values['TEST000001']['A'] = 999
    assert db.embedding('TEST000001')['TEST000001']['A'] == 0
    for ids in [[], ['TEST000001', 'TEST000001'], [1]]:
        with pytest.raises(ValueError):
            db.embedding(ids)
    with pytest.raises(KeyError):
        db.embedding('NONE000001')
    with pytest.raises(ValueError):
        db.search(' ')


def test_missing_values_are_inspectable_but_not_silently_imputed(tmp_path):
    db = database(tmp_path, record(values=['NA'] + list(range(1, 20))))
    entry = db['TEST000001']
    assert entry.missing_residues == ('A',) and not entry.complete
    assert db.search('synthetic', complete_only=True) == ()
    assert db.metadata('TEST000001')['entries'][0]['values']['A'] is None
    with pytest.raises(ValueError, match='TEST000001: A'):
        db.embedding('TEST000001')


@pytest.mark.parametrize('text,match', [
    ('', 'no records'),
    (record().replace('//\n', ''), 'unterminated'),
    (record()+record(), 'duplicate accession'),
    (record().replace('H TEST000001', 'H bad'), 'invalid AAindex1 accession'),
    (record().replace('R\n', 'H TEST000002\n'), 'duplicate H'),
    ('  orphan\n'+record(), 'continuation'),
    (record().replace('R\n', 'Z extra\n'), 'unrecognized'),
    (record().replace(HEADER, 'M rows = AR, cols = AR'), 'matrix records'),
    (record().replace(HEADER, 'I A/A R/K N/M D/F C/P Q/S E/T G/W H/Y I/V'), 'canonical'),
    (record(header='I A/L'), 'header'),
    (record(values=list(range(19))), 'expected 20'),
    (record(values=['nan']+list(range(1, 20))), 'nonfinite'),
    (record(values=['inf']+list(range(1, 20))), 'nonfinite'),
    (record(values=['oops']+list(range(1, 20))), 'invalid value'),
])
def test_malformed_records_fail_with_context(tmp_path, text, match):
    with pytest.raises(ValueError, match=match):
        database(tmp_path, text)


def test_digest_pin_rejects_changed_file(tmp_path):
    with pytest.raises(ValueError, match='SHA256 mismatch'):
        database(tmp_path, expected_sha256='0'*64)
    with pytest.raises(ValueError, match='64-character'):
        database(tmp_path, expected_sha256='bad')


def test_download_is_explicit_validates_before_writing_and_reuses_offline(tmp_path, monkeypatch):
    calls = []

    def retrieve(url, timeout):
        calls.append((url, timeout))
        return io.BytesIO(record().encode())

    monkeypatch.setattr('protspect.aaindex.urlopen', retrieve)
    path = tmp_path/'cache'/'aaindex1'
    db = download_aaindex1(path, timeout=3)
    assert calls == [(AAINDEX1_URL, 3)]
    assert db.source_url == AAINDEX1_URL
    assert db.sha256 == AAIndex.from_file(path).sha256
    with pytest.raises(FileExistsError):
        download_aaindex1(path)
    assert len(calls) == 1
    other = tmp_path/'invalid'
    with pytest.raises(ValueError, match='SHA256 mismatch'):
        download_aaindex1(other, expected_sha256='0'*64)
    assert not other.exists()
    calls.clear()
    with pytest.raises(ValueError, match='64-character'):
        download_aaindex1(other, expected_sha256='invalid')
    assert not calls
    monkeypatch.setattr('protspect.aaindex.urlopen', lambda *a, **k: io.BytesIO(b'bad'))
    with pytest.raises(ValueError):
        download_aaindex1(other)
    assert not other.exists()


def test_download_failure_preserves_existing_cache(tmp_path, monkeypatch):
    path = tmp_path/'data'
    def fail(*args, **kwargs):
        raise OSError('network unavailable')
    monkeypatch.setattr('protspect.aaindex.urlopen', fail)
    with pytest.raises(OSError, match='unavailable'):
        download_aaindex1(path)
    assert not path.exists()


def test_aaindex_channels_match_direct_dft_and_existing_feature_apis(tmp_path):
    db = database(tmp_path)
    embedding = db.embedding('TEST000001')
    item = SequenceFFT('ARLA')
    item.embed_sequence(embedding)
    item.compute_spectra(8)
    signal = np.array([0, 1, 10, 0, 0, 0, 0, 0])
    direct = np.exp(-2j*np.pi*np.outer(item.frequency, np.arange(8))) @ signal
    assert_allclose(item.spectra['TEST000001'], direct, atol=1e-12)
    dataset = SequenceSet(embedding=embedding, compute_cross_spectra=False)
    dataset.from_sequence_list(['ARLA'], spectrum_length=8)
    assert dataset.export_datatable().shape == (1, 4)
    paired = ComplexSet(embedding=embedding)
    paired.from_sequence_list(['AR'], ['LA'], spectrum_length=8)
    assert paired.export_datatable().shape == (1, 4)


def test_aaindex_estimator_cloning_pickling_and_snapshot(tmp_path):
    pytest.importorskip('sklearn')
    from sklearn.base import clone
    from protspect.sklearn import SpectralTransformer
    db = database(tmp_path)
    embedding = db.embedding('TEST000001')
    estimator = SpectralTransformer(embedding=embedding, spectrum_length=8)
    x = ['ARLA', 'LAAR']
    expected = estimator.fit_transform(x)
    assert_allclose(clone(estimator).fit_transform(x), expected)
    assert_allclose(pickle.loads(pickle.dumps(estimator)).transform(x), expected)
    embedding['TEST000001']['A'] = 999
    assert_allclose(estimator.transform(x), expected)
