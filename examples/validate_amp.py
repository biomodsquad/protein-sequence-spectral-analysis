"""Controlled ABP_GN correction audit using the previous classifier protocol.

Run from a checkout with the ``validation`` optional dependencies installed.
No encoder is trained. Only training data select classifier hyperparameters.
The supplied test set has been used in earlier research, so this is exploratory.
"""

import argparse
import hashlib
import json
import platform
import subprocess
import warnings
from collections import Counter
from pathlib import Path
from time import perf_counter

import numpy as np
import pandas as pd
import scipy
import sklearn
from scipy.fft import fft
from sklearn.base import clone
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (accuracy_score, average_precision_score, balanced_accuracy_score,
                             confusion_matrix, f1_score, matthews_corrcoef, roc_auc_score)
from sklearn.model_selection import train_test_split
from sklearn.neural_network import MLPClassifier
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from threadpoolctl import threadpool_limits

from protspect import SequenceFFT
from protspect.embeddings import b100


VARIANTS = {
    "legacy": (False, False),
    "fft_fix_only": (True, False),
    "matrix_fix_only": (False, True),
    "corrected": (True, True),
}


def load_data(path):
    """Read the same normalized Sequence/Class schema as the earlier benchmark."""
    frame = pd.read_csv(path)
    if not {"Sequence", "Class"}.issubset(frame):
        raise ValueError("CSV must contain Sequence and Class")
    if frame.Sequence.isna().any():
        raise ValueError("missing sequence")
    sequences = frame.Sequence.astype(str).str.upper().str.replace(r"\s+", "", regex=True).tolist()
    if any(not s or set(s) - set(b100) for s in sequences):
        raise ValueError("empty or noncanonical sequence")
    labels = frame.Class.map({"Negative": 0, "Positive": 1})
    if labels.isna().any() or labels.nunique() != 2:
        raise ValueError("both Positive and Negative classes are required")
    return sequences, labels.to_numpy(dtype=int)


def embedding_for(correct_matrix):
    """Return an independently copied original or corrected BLOSUM100 mapping."""
    embedding = {name: dict(values) for name, values in b100.items()}
    if not correct_matrix:
        embedding["P"]["A"] = 2
    return embedding


def feature_row(sequence, length, variant, cross):
    """Use package exports, replacing only bin order to reproduce legacy behavior."""
    correct_fft, correct_matrix = VARIANTS[variant]
    item = SequenceFFT(sequence)
    item.embed_sequence(embedding_for(correct_matrix))
    item.compute_spectra(length)
    if not correct_fft:
        # This is precisely the original unshifted FFT plus duplicated DC.
        item.frequency = np.linspace(-0.5, 0.5, length + 1)
        item.spectra = {}
        for name, vector in item.embeddings.items():
            values = fft(vector, n=length)
            item.spectra[name] = np.append(values, values[0])
    if cross:
        item.compute_cross_spectra()
    return item.compute_spectral_properties(("cross_spectra" if cross else "spectra",))


def metrics(labels, probabilities):
    """Compute the same thresholded and ranking metrics as previous experiments."""
    predicted = probabilities >= .5
    tn, fp, fn, tp = confusion_matrix(labels, predicted, labels=[0, 1]).ravel()
    return {
        "accuracy": float(accuracy_score(labels, predicted)),
        "balanced_accuracy": float(balanced_accuracy_score(labels, predicted)),
        "mcc": float(matthews_corrcoef(labels, predicted)),
        "f1": float(f1_score(labels, predicted)),
        "roc_auc": float(roc_auc_score(labels, probabilities)),
        "pr_auc": float(average_precision_score(labels, probabilities)),
        "sensitivity": float(tp / (tp + fn)), "specificity": float(tn / (tn + fp)),
    }


def candidates(kind, seed):
    """Reproduce the prior four-candidate LR/MLP search without encoder imports."""
    if kind == "linear":
        for c in (.01, .1, 1., 10.):
            yield {"C": c}, make_pipeline(StandardScaler(), LogisticRegression(
                C=c, class_weight="balanced", max_iter=3000, random_state=seed))
    else:
        for hidden, alpha in (((64,), 1e-4), ((64,), 1e-3),
                              ((128, 64), 1e-4), ((128, 64), 1e-3)):
            yield {"hidden": list(hidden), "alpha": alpha}, make_pipeline(
                StandardScaler(), MLPClassifier(hidden_layer_sizes=hidden, alpha=alpha,
                    activation="relu", early_stopping=True, validation_fraction=.15,
                    n_iter_no_change=20, max_iter=500, random_state=seed))


def fit_record(estimator, x, y):
    """Fit and preserve convergence warnings rather than hiding them."""
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", ConvergenceWarning)
        estimator.fit(x, y)
    return [str(w.message) for w in caught]


def select_classifier(kind, x, labels, test_x, test_labels, fit, validation, seed):
    """Select on validation BA/MCC/AP, refit on all training data, then test."""
    trials, best = [], None
    for settings, estimator in candidates(kind, seed):
        messages = fit_record(estimator, x[fit], labels[fit])
        score = metrics(labels[validation], estimator.predict_proba(x[validation])[:, 1])
        trials.append({"settings": settings, "validation": score, "warnings": messages})
        ranking = tuple(score[k] for k in ("balanced_accuracy", "mcc", "pr_auc"))
        if best is None or ranking > best[0]:
            best = ranking, settings, estimator, score
    _, settings, estimator, score = best
    final = clone(estimator)
    messages = fit_record(final, x, labels)
    probabilities = final.predict_proba(test_x)[:, 1]
    classifier = final.steps[-1][1]
    parameters = (sum(v.size for v in classifier.coefs_) +
                  sum(v.size for v in classifier.intercepts_)) if kind != "linear" else x.shape[1] + 1
    return {"selected_settings": settings, "validation": score,
            "test": metrics(test_labels, probabilities), "validation_trials": trials,
            "refit_warnings": messages, "classifier_parameters": parameters}, probabilities


def paired_interval(labels, old, new, seed, repeats=5000):
    """Stratified paired bootstrap CI for the change in balanced accuracy."""
    rng = np.random.default_rng(seed)
    delta = (new >= .5).astype(float) - (old >= .5).astype(float)
    positive, negative = delta[labels == 1], -delta[labels == 0]
    differences = np.empty(repeats)
    for i in range(repeats):
        differences[i] = .5 * (rng.choice(positive, len(positive)).mean() +
                               rng.choice(negative, len(negative)).mean())
    return {"delta": float(.5 * (positive.mean() + negative.mean())),
            "ci95": np.quantile(differences, [.025, .975]).tolist(),
            "bootstrap_repeats": repeats,
            "scope": "test examples only; excludes training/selection variability and homology"}


def digest(path):
    """Hash exact bytes for dataset and source provenance."""
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    """Generate audited features, classifiers, predictions, and provenance."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--train-csv", type=Path, required=True)
    parser.add_argument("--test-csv", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--fft-length", type=int, default=1024)
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument("--reference-report", type=Path)
    args = parser.parse_args()
    if args.fft_length < 2 or args.fft_length & (args.fft_length - 1):
        raise ValueError("use an even power-of-two length for the controlled legacy comparison")
    train, labels = load_data(args.train_csv)
    test, test_labels = load_data(args.test_csv)
    if set(train) & set(test):
        raise ValueError("train/test sequence overlap")
    if max(map(len, train + test)) > args.fft_length:
        raise ValueError("FFT length would truncate sequences")
    fit, validation = train_test_split(np.arange(len(train)), test_size=.2,
                                      random_state=args.seed, stratify=labels)
    output = args.output_dir
    output.mkdir(parents=True, exist_ok=False)
    np.savez_compressed(output / "split_indices.npz", fit=fit, validation=validation)
    report = {
        "task": "ABP_GN: Gram-negative antibacterial peptide classification",
        "protocol": {"seed": args.seed, "fft_length": args.fft_length,
            "selection": "stratified 80/20 training split; maximize BA, then MCC, then AP",
            "threshold": .5, "scaling": "StandardScaler fitted within each training partition",
            "normalization": "raw amplitudes; no length/energy normalization",
            "feature_types": ["spectra amplitude", "cross_spectra amplitude"],
            "selector": "frequency_select_physical: DC and positive dyadic frequencies",
            "test_status": "previously studied test set; exploratory correction audit",
            "uncertainty": "5000 paired stratified bootstrap replicates; no multiplicity correction"},
        "data": {}, "environment": {"python": platform.python_version(),
            "numpy": np.__version__, "scipy": scipy.__version__, "pandas": pd.__version__,
            "scikit_learn": sklearn.__version__, "blas_threads": 1},
        "source": {"script_sha256": digest(__file__), "package_files": {
            p.name: digest(p) for p in (Path(__file__).resolve().parents[1] / "protspect").glob('*.py')},
            "base_commit": subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
            "local_changes": True}, "models": {},
    }
    for name, path, sequences, y in (("train", args.train_csv, train, labels),
                                     ("test", args.test_csv, test, test_labels)):
        report["data"][name] = {"path": str(path.resolve()), "sha256": digest(path),
            "examples": len(y), "positive": int(y.sum()), "min_length": min(map(len, sequences)),
            "max_length": max(map(len, sequences)),
            "duplicate_sequences": sum(v - 1 for v in Counter(sequences).values())}
    predictions = pd.DataFrame({"Sequence": test, "Class": test_labels})
    all_sequences = train + test
    with threadpool_limits(limits=1):
        for cross in (False, True):
            kind_name = "cross_amplitude" if cross else "amplitude"
            for variant in VARIANTS:
                name = f"{kind_name}__{variant}"
                start = perf_counter()
                rows = [feature_row(sequence, args.fft_length, variant, cross)
                        for sequence in all_sequences]
                columns = list(rows[0])
                values = np.asarray([list(row.values()) for row in rows], dtype=np.float64)
                del rows
                if not np.isfinite(values).all():
                    raise ValueError("nonfinite spectral features")
                (output / f"{name}_columns.json").write_text(json.dumps(columns, indent=2) + '\n')
                np.savez_compressed(output / f"{name}_features.npz",
                                    train=values[:len(train)], test=values[len(train):])
                result = {"feature_count": len(columns), "feature_seconds": perf_counter() - start}
                for kind in ("linear", "nonlinear"):
                    start = perf_counter()
                    result[kind], probability = select_classifier(kind, values[:len(train)], labels,
                        values[len(train):], test_labels, fit, validation, args.seed)
                    result[kind]["fit_seconds"] = perf_counter() - start
                    predictions[f"{name}__{kind}"] = probability
                    print(name, kind, json.dumps(result[kind]["test"]), flush=True)
                report["models"][name] = result
                (output / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
                predictions.to_csv(output / 'test_predictions.csv', index=False)
    report['paired_comparisons'] = {}
    for feature in ('amplitude', 'cross_amplitude'):
        for kind in ('linear', 'nonlinear'):
            for variant in list(VARIANTS)[1:]:
                key = f'{feature}__{variant}__{kind}'
                report['paired_comparisons'][key] = paired_interval(test_labels,
                    predictions[f'{feature}__legacy__{kind}'].to_numpy(),
                    predictions[key].to_numpy(), args.seed)
    if args.reference_report:
        previous = json.loads(args.reference_report.read_text())
        report['historical_references'] = {
            'path': str(args.reference_report.resolve()), 'sha256': digest(args.reference_report),
            'status': 'previously saved scores; not retrained; different feature definitions',
            'models': {name: previous['models'][name] for name in (
                'blosum', 'blosum_dyadic_cross_spectral', 'esm2_8m')}}
    (output / 'report.json').write_text(json.dumps(report, indent=2) + '\n')


if __name__ == '__main__':
    main()
