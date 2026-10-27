# Reproducible CPP research example

Use [CPP_reproducible.ipynb](CPP_reproducible.ipynb), or run the same logic in
`cpp_workflow.py`. This replaces the historical exploration as the executable
research example. The original notebook in `FFTpro/examples/` is archived with
its source retained and all stale outputs removed. Its old personal import path,
obsolete arguments, unseeded choices, independently inferred grids, and saved
scores are not used by the corrected workflow.

## Isolated research environment

The tested environment is **CPython 3.12 on macOS**, with resolved versions in
[requirements-lock.txt](requirements-lock.txt). These are research dependencies;
none of MISTIC, notebooks, matplotlib, or scikit-learn is a core runtime
requirement. The `ml`/`validation` extras remain optional and separate.

MISTIC is the `mistic-svm` distribution from
[biomodsquad/interpretable-explainable-svms](https://github.com/biomodsquad/interpretable-explainable-svms),
version **0.2.0b1**, at immutable commit
[`8f115591921ed9b14fa03716246898694753780e`](https://github.com/biomodsquad/interpretable-explainable-svms/commit/8f115591921ed9b14fa03716246898694753780e).
The workflow verifies every installed MISTIC Python file against
[mistic-source.json](mistic-source.json), so a version string alone cannot hide
a dirty or different local checkout. The local validation used a pristine
`git archive` of that publicly verified commit, installed into the isolated
environment. It did not use the sibling working tree's uncommitted changes.

From the repository root:

```bash
python3.12 -m venv .venv-research
source .venv-research/bin/activate
python -m pip install setuptools==84.2.1 wheel==0.45.1
python -m pip install --no-build-isolation -r research/requirements-lock.txt
python -m pip install --no-build-isolation --no-deps -e .
python -m pip check
python -m ipykernel install --prefix .venv-research \
  --name protspect-research --display-name 'Protspect research (pinned)'
```

The lock captures the resolved research stack, not a claim of cross-platform
bitwise equivalence. A new Python/scientific-library stack needs revalidation.

## Prepare the complete website datasets

```bash
python research/prepare_cpp_data.py --output-dir research/data/mlcpp2
```

This downloads and verifies the complete official MLCPP 2.0 Layer1 snapshot,
then preserves every record, label, official FASTA identifier, and input order.
The committed source manifest pins the download hashes; changed files fail
verification. The raw and prepared data remain in an ignored local cache.
No historical CSV or subsetting decision is needed. See
[DATA_PROVENANCE.md](DATA_PROVENANCE.md) for source links and citations.

## Run headlessly

```bash
python research/cpp_workflow.py \
  --train-csv research/data/mlcpp2/training.csv \
  --test-csv research/data/mlcpp2/independent.csv \
  --output-dir validation/cpp_new_run
```

Or execute the notebook and save all outputs and figures:

```bash
python research/execute_notebook.py --output-dir validation/cpp_new_notebook_run
```

The output directory must be new; existing results are never silently reused
or overwritten. Jupyter needs permission to start a local kernel and bind local
ports. The executor sets workspace-local cache/configuration directories.
For an interactive notebook run, select the pinned kernel and choose a fresh
`OUTPUT` path in the setup cell. `--without-mistic` runs the baseline-only script
without importing MISTIC; it still needs the optional sklearn integration.

The research assets and datasets are not installed in the wheel. The source
checkout/distribution contains the preparation code and pinned source manifest;
obtain the actual data directly from the website. Both distribution formats
exclude raw downloads and prepared CSVs.

## Protocol fixed before test scoring

- Use seed 7, one BLAS thread, BLOSUM100 amplitude features, a predefined
  64-position FFT grid covering the training length range, the default dyadic/DC
  selector, and raw uncentered spectra. This gives 140 features. Test sequences
  exceeding this fixed grid raise an error; no truncation or grid resizing occurs.
- Reserve a stratified 20% **training-file** validation subset. Fit the
  SpectralTransformer and StandardScaler only on the other 80% during selection.
- Select the baseline RBF SVC from C = 1, 10 and gamma = 0.01, 0.1. Select the
  MISTIC ranking weight from 0, 0.5, 1 using the same external validation rows.
  Maximize balanced accuracy, then MCC, then average precision; ties keep the
  earlier candidate. No test score or explanation selects settings.
- MISTIC uses three independently partitioned member groups, 20% internal
  validation, balanced precomputed-kernel SVCs, shared features, training-set
  ranking, 50% backward removal, initial tuning followed by scaled-gamma search,
  and final retuning/knee selection. Its internal splits reuse the scaler fitted
  on their outer training input; their scores are **not fold-pure or unbiased
  validation estimates**. The outer validation and final test are excluded from
  every such fitting/search input.
- Refit the chosen algorithm and preprocessing on all training records, then
  evaluate all official test records except four exact training duplicates
  (2,337 retained: 153 positive, 2,184 negative).
  MISTIC explicitly uses member-ensemble decisions and its training-derived
  cutoff, not the newer default unified predictor or test-tuned thresholds.
- Score both predeclared examples. Do not select or tune a winner using test
  results. This compact search has a different budget/protocol from the original
  exploratory notebook and is not a reproduction of its numerical results.

ROC-AUC and average precision use decision scores, not calibrated probabilities.
Brier/calibration curve fields are missing (`null`) because `probability=False`;
these are not zero-valued performance measurements. The report captures warnings.

## Explanations and saved evidence

The workflow explains 32 deterministic test examples, sampled without labels,
only after fitting. Integrated gradients use a zero standardized-feature vector
as baseline and explain the mean member-ensemble decision. The output concerns
spectral features, not residues or a causal biological mechanism.

Each run stores official source URLs/raw and prepared-file hashes, versions,
official FASTA identifiers, original-row split IDs, exact
exclusions, selected settings, all external-validation trials, selected features,
training-derived cutoff, feature-removal curve, exact ordered feature metadata,
test decisions/predictions, explanation sample IDs, and attribution summaries.
Raw peptide sequences are not duplicated into the reports. Original source data
and the pinned stack are required for refitting; a model pickle is not shipped.

## Check the example

```bash
python -m pytest -q research/tests
```

The protocol tests cover overlap/label handling, FASTA comparison, source row
preservation, strict JSON missing statistics, deterministic baseline replay, and
refusal to mix runs. They supplement the package tests in `tests/`.
The optional manual GitHub Actions workflow `research.yml` installs the pinned
research environment, prepares the complete official website datasets, runs
these checks, and executes the notebook. It has not
been dispatched as part of this local review.
