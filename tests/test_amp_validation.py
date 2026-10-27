"""Check that the correction audit reproduces the actual upstream algorithm."""

import importlib.util
import subprocess
from pathlib import Path

import numpy as np
import pytest
from numpy.testing import assert_allclose

pytest.importorskip("sklearn")

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('validate_amp', ROOT / 'examples/validate_amp.py')
audit = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(audit)


def test_legacy_reproduction_against_original_source():
    if not (ROOT / ".git").exists():
        pytest.skip("original-source comparison requires a git checkout")
    if subprocess.run(["git", "cat-file", "-e",
            "b31d4da5fbea0fdcc3c948429ddade0a7dfba3d2"], cwd=ROOT,
            capture_output=True).returncode:
        pytest.skip("original commit is unavailable in this checkout")
    source = subprocess.check_output(
        ['git', 'show', 'b31d4da5fbea0fdcc3c948429ddade0a7dfba3d2:FFTpro/features.py'],
        cwd=ROOT, text=True).replace('from FFTpro.', 'from protspect.')
    namespace = {}
    exec(compile(source, '<upstream-features>', 'exec'), namespace)
    for sequence in ['ACDE', 'PPAAKRW', 'ACDEFGHIKLMNPQRSTVWY']:
        item = namespace['SequenceFFT'](sequence)
        item.embed_sequence(audit.embedding_for(False))
        item.compute_spectra(1024)
        item.compute_cross_spectra()
        for cross in (False, True):
            original = item.compute_spectral_properties(
                ['cross_spectra' if cross else 'spectra'])
            actual = audit.feature_row(sequence, 1024, 'legacy', cross)
            assert list(original) == list(actual)
            assert_allclose(list(original.values()), list(actual.values()), atol=1e-10)


def test_matrix_fix_has_no_effect_without_alanine():
    for cross in (False, True):
        old = audit.feature_row('KRWPP', 1024, 'legacy', cross)
        new = audit.feature_row('KRWPP', 1024, 'matrix_fix_only', cross)
        assert_allclose(list(old.values()), list(new.values()))


def test_bootstrap_identical_predictions_has_zero_delta():
    result = audit.paired_interval(np.array([0, 1, 0, 1]),
        np.array([.2, .8, .7, .6]), np.array([.2, .8, .7, .6]), 7, repeats=20)
    assert result['delta'] == 0
    assert result['ci95'] == [0, 0]


def test_missing_sequence_rejected(tmp_path):
    path = tmp_path / 'bad.csv'
    path.write_text('Sequence,Class\n,Positive\nAC,Negative\n')
    with pytest.raises(ValueError, match='missing'):
        audit.load_data(path)
