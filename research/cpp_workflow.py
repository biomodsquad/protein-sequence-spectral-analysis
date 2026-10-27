"""Reproducible CPP research example with an untouched outer validation subset.

This is a corrected, compact replacement for the historical exploration,
not an attempt to recover its old numerical results or a new independent benchmark.
"""

import argparse
from contextlib import redirect_stdout
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform
import warnings

import numpy as np
import pandas as pd
from sklearn.base import clone
from sklearn.metrics import (accuracy_score, average_precision_score, balanced_accuracy_score,
                             confusion_matrix, f1_score, matthews_corrcoef, roc_auc_score)
from sklearn.model_selection import train_test_split
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from threadpoolctl import threadpool_limits

from protspect.sklearn import SpectralTransformer

ROOT = Path(__file__).resolve().parents[1]
SEED = 7
FFT_LENGTH = 64
C_VALUES = (1., 10.)
GAMMA_VALUES = (.01, .1)
RANK_WEIGHTS = (0., .5, 1.)


def sha256(path):
    """Return a SHA-256 digest of exact source bytes."""
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_table(path):
    """Load canonical, uniquely labeled sequences without silently imputing data."""
    table = pd.read_csv(path)
    if not {'Sequence', 'Class'}.issubset(table.columns) or table[['Sequence', 'Class']].isna().any().any():
        raise ValueError('CSV requires nonmissing Sequence and Class columns')
    if 'SourceID' in table:
        if table.SourceID.isna().any() or table.SourceID.duplicated().any():
            raise ValueError('SourceID values must be nonmissing and unique')
        table = table[['Sequence', 'Class', 'SourceID']].copy()
        table = table.rename(columns={'SourceID': 'source_id'})
    else:
        table = table[['Sequence', 'Class']].copy()
    table['source_row'] = np.arange(len(table))
    table.Sequence = table.Sequence.astype(str).str.upper().str.replace(r'\s+', '', regex=True)
    if any(not sequence or set(sequence) - set('ACDEFGHIKLMNPQRSTVWY') for sequence in table.Sequence):
        raise ValueError('CSV contains empty or noncanonical sequences')
    labels = table.Class.map({'Negative': 0, 'Positive': 1})
    if labels.isna().any() or labels.nunique() != 2:
        raise ValueError('CSV must contain both Positive and Negative classes')
    if table.Sequence.duplicated().any():
        raise ValueError('duplicate sequences require explicit resolution before analysis')
    table.Class = labels.astype(int)
    return table


def clean_test_overlap(train, test):
    """Exclude exact test duplicates of training sequences and reject label conflicts."""
    labels = train.set_index('Sequence').Class.to_dict()
    overlap = test.Sequence.isin(labels)
    excluded = test.loc[overlap].copy()
    for row in excluded.itertuples():
        if labels[row.Sequence] != row.Class:
            raise ValueError('train/test duplicates have conflicting class labels')
    clean = test.loc[~overlap].copy()
    if clean.Class.nunique() != 2:
        raise ValueError('cleaned test set must contain both classes')
    audit = [{'test_source_row': int(row.source_row), 'class': int(row.Class),
              'sequence_sha256': hashlib.sha256(row.Sequence.encode()).hexdigest(),
              **({'source_id': row.source_id} if hasattr(row, 'source_id') else {})}
             for row in excluded.itertuples()]
    return clean.reset_index(drop=True), audit


def metrics(labels, decision, predicted):
    """Report class-balanced and ranking metrics, with AP from decision scores."""
    tn, fp, fn, tp = confusion_matrix(labels, predicted, labels=[0, 1]).ravel()
    return {'accuracy': float(accuracy_score(labels, predicted)),
            'balanced_accuracy': float(balanced_accuracy_score(labels, predicted)),
            'mcc': float(matthews_corrcoef(labels, predicted)),
            'f1': float(f1_score(labels, predicted, zero_division=0)),
            'roc_auc': float(roc_auc_score(labels, decision)),
            'average_precision': float(average_precision_score(labels, decision)),
            'sensitivity': float(tp / (tp + fn)), 'specificity': float(tn / (tn + fp)),
            'confusion_matrix': [[int(tn), int(fp)], [int(fn), int(tp)]]}


def ranking(result):
    """Choose candidates using outer validation alone, with deterministic tie breaks."""
    return tuple(result[key] for key in ('balanced_accuracy', 'mcc', 'average_precision'))


def verify_mistic():
    """Reject a dirty or different MISTIC implementation despite a matching version string."""
    import mistic
    expected = json.loads((ROOT / 'research/mistic-source.json').read_text())
    source = Path(mistic.__file__).resolve().parent
    actual = {str(path.relative_to(source)): sha256(path) for path in sorted(source.rglob('*.py'))}
    if actual != expected['files'] or mistic.__version__ != expected['version']:
        raise RuntimeError('MISTIC source differs from pinned commit; install research/requirements-lock.txt')
    return expected


def mistic_fit(values, labels, weight, log_path):
    """Fit a deterministic training-only feature search using a compact parameter grid."""
    from mistic import combined_rank, cvSet, kernelWrapper, paramSet, score_svc, svmSet
    splits = cvSet(values, labels, num_feature_medoids=0)
    splits.independent(num_sets=3, validation_size=.2, random_seed=SEED)
    model = svmSet(SVC(kernel='precomputed', class_weight='balanced', probability=False,
                       random_state=SEED), splits, score_svc().score,
                   kernel=kernelWrapper('rbf'), feature_set_policy='shared',
                   kernel_feature_selection='shared')
    grid = [paramSet({'C': c}, {'gamma': gamma}) for c in C_VALUES for gamma in GAMMA_VALUES]
    with log_path.open('w') as stream, redirect_stdout(stream):
        model.greedy_backward_selection(grid, reduction_factor=.5,
            feature_ranker=combined_rank(weight=weight, random_seed=SEED).compute,
            set_for_rank='train', tune_models_each_step=False, post_find_knee=True)
    return model


def mistic_output(model, values):
    """Use explicit member-ensemble scores and the training-derived cutoff consistently."""
    decision = model.decision_function(values, prediction_mode='set')
    predicted = model.predict(values, prediction_mode='set')
    return np.asarray(decision), np.asarray(predicted, dtype=int)


def jsonable(value):
    """Convert library result containers to stable JSON values."""
    if isinstance(value, dict):
        return {str(key): jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [jsonable(item) for item in value]
    if isinstance(value, np.generic):
        return jsonable(value.item())
    if isinstance(value, float) and not np.isfinite(value):
        # Unavailable library statistics (e.g. probability metrics with probability=False) are missing.
        return None
    return value


def run(train_csv, test_csv, output_dir, with_mistic=True):
    """Run the complete analysis and save source/data hashes, splits, trials, and predictions.

    Select settings on a stratified 20% holdout of the training file, refit the
    selected algorithms on all training data, then evaluate the exact-disjoint
    test subset. The source test file is never changed. No homology claim is made.
    """
    train_csv, test_csv, output = Path(train_csv), Path(test_csv), Path(output_dir)
    source_train, source_test = read_table(train_csv), read_table(test_csv)
    test, exclusions = clean_test_overlap(source_train, source_test)
    if max(source_train.Sequence.str.len().max(), test.Sequence.str.len().max()) > FFT_LENGTH:
        raise ValueError('the predefined 64-position CPP grid would truncate a sequence')
    mistic_source = verify_mistic() if with_mistic else None
    output.mkdir(parents=True, exist_ok=False)
    np.random.seed(SEED)
    fit, validation = train_test_split(np.arange(len(source_train)), test_size=.2,
        stratify=source_train.Class.to_numpy(), random_state=SEED)
    split_rows = pd.DataFrame({'source_row': source_train.source_row, 'partition': 'fit'})
    split_rows.loc[validation, 'partition'] = 'validation'
    if 'source_id' in source_train:
        split_rows['source_id'] = source_train.source_id
    split_rows.to_csv(output/'training_split.csv', index=False)
    pd.DataFrame(exclusions).to_csv(output/'excluded_test_rows.csv', index=False)
    test[[name for name in ('source_row', 'source_id') if name in test]].to_csv(
        output/'retained_test_rows.csv', index=False)
    y = source_train.Class.to_numpy()
    y_test = test.Class.to_numpy()
    preprocessing = make_pipeline(SpectralTransformer(spectrum_length=FFT_LENGTH), StandardScaler())
    report = {
        'purpose': 'compact reproducible research example; not historical result reproduction',
        'protocol': {'seed': SEED, 'fft_length': FFT_LENGTH, 'outer_validation_fraction': .2,
            'selection': 'outer validation balanced accuracy, MCC, then average precision',
            'grid_C': list(C_VALUES), 'grid_gamma': list(GAMMA_VALUES),
            'rank_weights': list(RANK_WEIGHTS), 'blas_threads': 1,
            'preprocessing': 'raw BLOSUM100 amplitudes, default dyadic/DC selector; fit scaler only on fit rows',
            'mistic_internal_validation': 'training-only heuristic for tuning and feature removal; its scores are not unbiased performance estimates',
            'mistic_search': '3 independent member groups; 20% internal validation; 50% removal per step; scaled-gamma search then retune; knee enabled',
            'mistic_prediction': 'member mean decision, training-derived decision cutoff, no voting',
            'test_usage': 'evaluate only after settings fixed; exclude exact train/test duplicates',
            'limits': 'previously explored data; no family-disjoint split; no claim of independent biological validation'},
        'data': {'train': {'file': train_csv.name, 'sha256': sha256(train_csv),
                          'rows': len(source_train), 'positive': int(y.sum())},
                 'test_original': {'file': test_csv.name, 'sha256': sha256(test_csv),
                          'rows': len(source_test), 'positive': int(source_test.Class.sum())},
                 'test_exact_disjoint': {'rows': len(test), 'positive': int(y_test.sum()),
                                        'excluded_exact_duplicates': exclusions}},
        'environment': {'python': platform.python_version(), 'system': platform.system(),
            'packages': {name: importlib.metadata.version(name) for name in (
                'numpy', 'scipy', 'pandas', 'scikit-learn', 'matplotlib', 'threadpoolctl', 'joblib')}},
        'source': {'workflow_sha256': sha256(__file__),
                   'protspect_files': {p.name: sha256(p) for p in sorted((ROOT/'protspect').glob('*.py'))},
                   'mistic': mistic_source}, 'models': {},
    }
    predictions = pd.DataFrame({'source_row': test.source_row, 'Class': y_test})
    if 'source_id' in test:
        predictions['source_id'] = test.source_id
    manifest_path = train_csv.parent / 'source.json'
    if manifest_path.exists() and train_csv.parent == test_csv.parent:
        provenance = json.loads(manifest_path.read_text())
        if (sha256(train_csv) != provenance['datasets']['training']['csv_sha256'] or
                sha256(test_csv) != provenance['datasets']['independent']['csv_sha256']):
            raise ValueError('prepared data no longer match their official source manifest')
        report['data']['official_source'] = provenance
    with warnings.catch_warnings(record=True) as caught, threadpool_limits(limits=1):
        warnings.simplefilter('always')
        fit_x = preprocessing.fit_transform(source_train.Sequence.iloc[fit].tolist())
        validation_x = preprocessing.transform(source_train.Sequence.iloc[validation].tolist())
        # Neither test labels nor test features enter any candidate selection.
        svc_trials, best_svc = [], None
        for c in C_VALUES:
            for gamma in GAMMA_VALUES:
                model = SVC(C=c, gamma=gamma, kernel='rbf', class_weight='balanced', random_state=SEED)
                model.fit(fit_x, y[fit])
                result = metrics(y[validation], model.decision_function(validation_x), model.predict(validation_x))
                settings = {'C': c, 'gamma': gamma}
                svc_trials.append({'settings': settings, 'validation': result})
                if best_svc is None or ranking(result) > best_svc[0]:
                    best_svc = ranking(result), settings, result
        # The full-training scaler and schema are refitted only after selection.
        final_preprocessing = clone(preprocessing)
        all_train_x = final_preprocessing.fit_transform(source_train.Sequence.tolist())
        (output/'feature_metadata.json').write_text(json.dumps(
            final_preprocessing.steps[0][1].get_feature_metadata(), indent=2) + '\n')
        test_x = final_preprocessing.transform(test.Sequence.tolist())
        final_svc = SVC(**best_svc[1], kernel='rbf', class_weight='balanced', random_state=SEED)
        final_svc.fit(all_train_x, y)
        svc_decision, svc_prediction = final_svc.decision_function(test_x), final_svc.predict(test_x)
        report['models']['rbf_svc'] = {'selected_settings': best_svc[1],
            'validation': best_svc[2], 'validation_trials': svc_trials,
            'test': metrics(y_test, svc_decision, svc_prediction)}
        predictions['rbf_svc_decision'] = svc_decision
        predictions['rbf_svc_prediction'] = svc_prediction
        if with_mistic:
            trials, best = [], None
            for weight in RANK_WEIGHTS:
                model = mistic_fit(fit_x, y[fit], weight, output/f'mistic_weight_{weight}.log')
                decision, predicted = mistic_output(model, validation_x)
                result = metrics(y[validation], decision, predicted)
                trials.append({'weight': weight, 'validation': result,
                               'selected_feature_count': len(model.features)})
                if best is None or ranking(result) > best[0]:
                    best = ranking(result), weight, result
            final_mistic = mistic_fit(all_train_x, y, best[1], output/'mistic_final.log')
            decision, predicted = mistic_output(final_mistic, test_x)
            columns = final_preprocessing.steps[0][1].get_feature_names_out()
            selected = np.asarray(final_mistic.features, dtype=int)
            report['models']['mistic'] = {'selected_rank_weight': best[1],
                'validation': best[2], 'validation_trials': trials,
                'test': metrics(y_test, decision, predicted),
                'selected_features': columns[selected].tolist(),
                'selected_feature_indices': selected.tolist(),
                'parameters': {'model': final_mistic.parameters_.model, 'kernel': final_mistic.parameters_.kernel},
                'decision_cutoff': float(final_mistic.decision_value_cutoff_),
                'selection_curve': jsonable(final_mistic.feature_performance_),
                'curve_missing_values': 'unavailable statistics, including disabled probability metrics, are serialized as null'}
            predictions['mistic_decision'], predictions['mistic_prediction'] = decision, predicted
            # Explanation rows are chosen without using test labels; no choices are
            # fed back into model selection. Attributions describe spectral features.
            rows = np.sort(np.random.default_rng(SEED).choice(len(test), min(32, len(test)), replace=False))
            attributions = final_mistic.integrated_gradient(test_x[rows],
                reference_point=np.zeros(all_train_x.shape[1]), num_steps=32)
            if attributions.shape[1] != len(selected):
                raise ValueError('unexpected MISTIC attribution column schema')
            pd.DataFrame({'feature': columns[selected],
                          'mean_absolute_integrated_gradient': np.abs(attributions).mean(axis=0)}
                         ).sort_values('mean_absolute_integrated_gradient', ascending=False).to_csv(
                             output/'attribution_summary.csv', index=False)
            pd.DataFrame({'source_row': test.source_row.iloc[rows]}).to_csv(
                output/'attribution_rows.csv', index=False)
    report['warnings'] = sorted(set(f'{type(item.message).__name__}: {item.message}' for item in caught))
    predictions.to_csv(output/'test_predictions.csv', index=False, float_format='%.17g')
    (output/'report.json').write_text(json.dumps(jsonable(report), indent=2, allow_nan=False) + '\n')
    return report


def main():
    """Run the compact experiment without notebook-specific globals or paths."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--train-csv', type=Path, required=True)
    parser.add_argument('--test-csv', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--without-mistic', action='store_true')
    args = parser.parse_args()
    result = run(args.train_csv, args.test_csv, args.output_dir, not args.without_mistic)
    print(json.dumps({name: model['test'] for name, model in result['models'].items()}, indent=2))


if __name__ == '__main__':
    main()
