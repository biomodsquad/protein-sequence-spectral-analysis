# ABP_GN four-transform MLP baseline

Executed 2026-10-27 using [the editable notebook](../research/AMP_transform_MLP.ipynb).
Raw data and row-indexed predictions stay local. This test set was previously
studied; these are exploratory single-seed results, not evidence of significant
transform differences.

| Transform | Validation BA | Test BA | Test MCC | Test ROC-AUC | Test average precision |
| --- | ---: | ---: | ---: | ---: | ---: |
| fourier | 0.8793 | 76.44% | 0.529 | 0.834 | 0.749 |
| z | 0.8676 | 76.38% | 0.525 | 0.845 | 0.782 |
| mellin | 0.7957 | 74.66% | 0.488 | 0.835 | 0.787 |
| wavelet | 0.8517 | 74.72% | 0.493 | 0.818 | 0.719 |

All four use corrected BLOSUM100, 20 channels × 256 coefficient magnitudes,
5,120 features, and 327,809 MLP parameters. No centering or per-sequence
normalization. Z radius 1.01; finite Mellin c=0; Haar wavelet level 4 with
periodization. Wavelet magnitudes discard signs. Fourier uses every shifted bin,
including DC and redundant conjugate magnitudes; no dyadic selector or cross
spectra. These results therefore do not reproduce the earlier FFT correction audit.

Preserved 3,562 supplied training and 1,973 test peptides. Exact file hashes and
seed-7 stratified 80/20 training validation indices match the earlier audit.
One common candidate: ReLU hidden layer (64,), alpha=0.001, Adam learning rate
0.001, batch size 128, maximum 300 epochs, training-only internal early stopping
(fraction 0.15, patience 15), seed 7. StandardScaler fits only each fitting
partition. Validation selects configurations before full-training refit; the test
set never selects hyperparameters or thresholds. Threshold fixed at 0.5.

Runtime approximately 2.1 minutes with one numerical thread. No validation or
refit convergence warnings. Mellin extraction took approximately 90 seconds;
other transforms each took 6–7 seconds. All notebook code cells executed, metrics
and predictions were checked, both figures visually inspected, and the
23 transform tests passed. Defaults were not adjusted after inspecting test scores.

See [metrics](amp_transform_mlp_baseline/metrics.csv) and
[controls, source/data hashes, trials, and runtime versions](amp_transform_mlp_baseline/report.json).
Run the notebook to regenerate full feature names, predictions, curves, and
split indices in ignored `validation/amp_transform_mlp/`.

Further comparison should include multiple model seeds and a family-disjoint or
new external test set. This run does not quantify homology or dataset uncertainty.
