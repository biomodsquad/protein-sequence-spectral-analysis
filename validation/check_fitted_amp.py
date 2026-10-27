"""Check fitted preprocessing on the existing ABP_GN files without tuning on test labels."""

import argparse
import hashlib
import json
from pathlib import Path
from time import perf_counter

import numpy as np
import pandas as pd

from protspect.sklearn import SpectralTransformer


def main():
    """Fit on the supplied training sequences and verify immutable test transforms."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--train-csv', type=Path, required=True)
    parser.add_argument('--test-csv', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    output = args.output_dir
    output.mkdir(exist_ok=True)
    paths = {'train': args.train_csv, 'test': args.test_csv}
    sequences = {name: pd.read_csv(path).Sequence.str.upper().str.replace(
        r'\s+', '', regex=True).tolist() for name, path in paths.items()}
    start = perf_counter()
    transformer = SpectralTransformer(spectrum_length=1024, normalization='length',
                                      standardize_channels=True)
    train = transformer.fit_transform(sequences['train'])
    metadata = transformer.get_feature_metadata()
    test = transformer.transform(sequences['test'])
    assert transformer.get_feature_metadata() == metadata
    np.testing.assert_allclose(transformer.transform(sequences['test'][:1]), test[:1])
    assert train.shape == (len(sequences['train']), 220)
    assert test.shape == (len(sequences['test']), 220)
    assert np.isfinite(train).all() and np.isfinite(test).all()
    (output / 'feature_metadata.json').write_text(json.dumps(metadata, indent=2) + '\n')
    report = {
        'purpose': 'integration check only; no classifier or normalization-policy selection',
        'preprocessing_fit': 'training sequences only; test labels not accessed',
        'train_shape': list(train.shape), 'test_shape': list(test.shape),
        'test_transform_preserves_fitted_metadata': True,
        'single_sequence_and_batch_outputs_agree': True,
        'finite_features': True, 'seconds': perf_counter() - start,
        'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'data_sha256': {name: hashlib.sha256(path.read_bytes()).hexdigest()
                        for name, path in paths.items()},
    }
    (output / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
