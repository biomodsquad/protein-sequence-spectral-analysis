# Scientific validation

The ABP_GN correction audit compares original and corrected sequence-feature
calculations on the same Gram-negative antibacterial peptide dataset and
classifier protocol previously used for the Fourier, wavelet, and ESM models.

From the repository root:

```bash
python -m pip install -e '.[dev,validation]'
python examples/validate_amp.py \
  --train-csv /path/to/ABP_GN_train.csv \
  --test-csv /path/to/ABP_GN_test.csv \
  --output-dir validation/new_abp_gn_run \
  --reference-report /path/to/previous/report.json
```

The output directory must not already exist. Defaults are FFT length 1024,
seed 7, a stratified 20% validation partition of the supplied training set,
and the prior four-candidate searches for linear and nonlinear classifiers.
The optional reference report copies historical scores for context; no pretrained
encoders or external checkpoints are required for this audit.

The 2×2 correction design tests neither fix, FFT-bin alignment only, the
BLOSUM100 P-to-A fix only, and both fixes, for two feature types: sequence
amplitudes and upper-triangle channel cross-spectrum amplitudes. These are the
package's original/default magnitude exports. They differ from the real/imaginary
cross features used in the strongest earlier benchmark. Raw amplitudes remain
unnormalized. The shared grid and feature count are fixed across variants.

Each run records dataset SHA-256 hashes, split indices, environment versions,
feature/source hashes, feature column names, raw feature matrices, hyperparameter
trials, convergence warnings, selected settings, predictions, and metrics.
Feature matrices are excluded from Git because they are reproducible and large.
Source data are read from their existing location rather than redistributed.
The raw data's upstream citation and redistribution terms still need confirmation.

The paired 95% bootstrap intervals describe changes in test balanced accuracy,
resampling positives and negatives separately and preserving paired predictions.
They do not measure variability from model fitting, hyperparameter selection,
or homologous peptide families; multiple contrasts are unadjusted. The test set
has already been inspected in prior research, and is not a new independent
holdout. No sequence-identity or family-disjoint evaluation was performed.

## Fitted preprocessing integration check

`check_fitted_amp.py` fits the new transformer on training sequences only,
checks finite train/test matrices and invariant metadata, and verifies that a
sequence has the same output alone or in a batch. It does not train a classifier
or select normalization based on test outcomes:

```bash
python validation/check_fitted_amp.py \
  --train-csv /path/to/ABP_GN_train.csv \
  --test-csv /path/to/ABP_GN_test.csv \
  --output-dir validation/normalization_amp_smoke
```

The report and exact feature provenance are saved in the output directory.
See `docs/PREPROCESSING.md` for policy definitions and pipeline usage.
