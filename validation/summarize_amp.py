"""Create human-readable tables and a scientific figure from the completed audit."""

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.switch_backend('Agg')

ROOT = Path(__file__).resolve().parent
RUN = ROOT / 'abp_gn_fft1024_seed7'
r = json.loads((RUN / 'report.json').read_text())
rows = []
for name, result in r['models'].items():
    for classifier in ('linear', 'nonlinear'):
        rows.append({'representation': name, 'classifier': classifier,
                     'features': result['feature_count'],
                     'classifier_parameters': result[classifier]['classifier_parameters'],
                     'selected_settings': json.dumps(result[classifier]['selected_settings']),
                     **result[classifier]['test']})
pd.DataFrame(rows).to_csv(RUN / 'metrics.csv', index=False)

variants = ('legacy', 'fft_fix_only', 'matrix_fix_only', 'corrected')
labels = ('Original', 'FFT fix only', 'Matrix fix only', 'Both fixes')
fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True)
for ax, feature, title in zip(axes, ('amplitude', 'cross_amplitude'),
                             ('Sequence amplitudes · 220 features',
                              'Cross-spectrum amplitudes · 2,310 features')):
    for offset, classifier, marker in ((-.12, 'linear', 'o'), (.12, 'nonlinear', 's')):
        values = [100*r['models'][f'{feature}__{v}'][classifier]['test']['balanced_accuracy']
                  for v in variants]
        ax.plot(np.arange(4) + offset, values, marker=marker, linewidth=1.3,
                label='Logistic regression' if classifier == 'linear' else 'MLP')
    ax.set_xticks(np.arange(4), labels, rotation=18)
    ax.set_ylim(67, 78)
    ax.set_title(title, fontsize=11)
    ax.set_ylabel('Test balanced accuracy (%)')
    ax.grid(axis='y', alpha=.2)
    ax.spines[['top', 'right']].set_visible(False)
axes[0].legend(frameon=False, fontsize=9)
fig.suptitle('ABP_GN correction audit · fixed split and validation protocol', fontsize=13)
fig.savefig(RUN/'correction_audit.png', dpi=180)
plt.close(fig)

lines = ['# ABP_GN: scientific validation of the corrections', '',
    'Run date: 2026-10-27. This is an exploratory audit on the previously studied split.', '',
    '## Main results', '',
    'Classifier hyperparameters were selected separately for each representation using only',
    'the stratified training validation partition. Results use the full training refit',
    'and the supplied test set at probability threshold 0.5.', '',
    '| Features | Classifier | Original BA | Corrected BA | Change (pp), paired 95% CI |',
    '| --- | --- | ---: | ---: | ---: |']
for feature in ('amplitude', 'cross_amplitude'):
    for classifier in ('linear', 'nonlinear'):
        old=r['models'][f'{feature}__legacy'][classifier]['test']['balanced_accuracy']
        new=r['models'][f'{feature}__corrected'][classifier]['test']['balanced_accuracy']
        interval=r['paired_comparisons'][f'{feature}__corrected__{classifier}']
        lo,hi=np.array(interval['ci95'])*100
        lines.append(f'| {feature} | {classifier} | {100*old:.2f}% | {100*new:.2f}% | '
                     f'{100*(new-old):+.2f} [{lo:+.2f}, {hi:+.2f}] |')
lines += ['', '## Nonlinear classifier ablation', '',
          '| Features | Variant | BA | MCC | ROC-AUC | Average precision | Hidden units; alpha |',
          '| --- | --- | ---: | ---: | ---: | ---: | --- |']
for feature in ('amplitude','cross_amplitude'):
    for variant in variants:
        model=r['models'][f'{feature}__{variant}']['nonlinear']
        m=model['test']
        s=model['selected_settings']
        lines.append(f'| {feature} | {variant} | {100*m["balanced_accuracy"]:.2f}% | '
                     f'{m["mcc"]:.3f} | {m["roc_auc"]:.3f} | {m["pr_auc"]:.3f} | '
                     f'{s["hidden"]}; {s["alpha"]} |')
lines += ['', '## Historical reference scores', '',
          'These are saved results, not newly retrained models. Test sequence order and labels',
          'were checked against the prior predictions and match. Feature definitions differ.', '',
          '| Earlier representation | Classifier | BA | MCC | ROC-AUC | Average precision |',
          '| --- | --- | ---: | ---: | ---: | ---: |']
for name in ('blosum', 'blosum_dyadic_cross_spectral', 'esm2_8m'):
    for classifier in ('linear', 'nonlinear'):
        m=r['historical_references']['models'][name][classifier]['test']
        lines.append(f'| {name} | {classifier} | {100*m["balanced_accuracy"]:.2f}% | '
                     f'{m["mcc"]:.3f} | {m["roc_auc"]:.3f} | {m["pr_auc"]:.3f} |')
lines += ['', '## Interpretation and limits', '',
    '- With the nonlinear classifier, corrected sequence amplitudes score 74.99%',
    '  balanced accuracy versus 74.22% originally; corrected cross amplitudes score',
    '  74.93% versus 75.30%. Both paired intervals include zero, so this run does not',
    '  establish a nonlinear performance gain or loss. Linear cross amplitudes improve',
    '  by 2.39 percentage points (unadjusted paired 95% interval +0.31 to +4.48).',
    '- No classifier convergence warnings were recorded.',
    '- The corrections restore the stated mathematical frequency semantics and the',
    '  reference substitution scores. Their validity does not depend on improving prediction.',
    '- Corrected sequence amplitudes and corrected cross-spectrum amplitudes must not',
    '  be treated as the earlier real/imaginary cross-spectrum representation.',
    '  In particular, `abs(X_i * conjugate(X_j)) = abs(X_i) * abs(X_j)`: the',
    '  amplitude-only cross representation contains quadratic magnitude interactions,',
    '  but loses relative phase. The earlier strongest cross-spectrum model retains',
    '  real and imaginary components, excludes DC, and uses a different frequency set.',
    '- The original bug changed which coefficients were exported, not the FFT itself.',
    '  Existing classifiers can exploit those coefficients even when their frequency labels',
    '  are wrong. A loss of predictive performance does not justify retaining wrong labels.',
    '- Both changes affect feature values, so recomputation and retraining remain required.',
    '- Data: 3,562 training peptides (1,930 positive), 1,973 test peptides (826 positive).',
    '  Maximum length is 255 residues; FFT length is fixed at 1,024, with no truncation.',
    '  No exact sequence duplicates occur within or across the supplied splits.',
    '- Seed 7, 713 validation examples, four hyperparameter candidates per classifier,',
    '  training-only StandardScaler, raw amplitudes, and the default dyadic/DC selector.',
    '  Feature counts remain identical across the four variants.',
    '- The paired bootstrap intervals use 5,000 stratified test resamples and do not',
    '  include fitting variability, model-selection variability, or homology structure.',
    '  No correction for multiple contrasts was applied. The test set was already used',
    '  in previous experiments. Treat small differences as tentative.',
    '- Next confirmation should use repeated training seeds and a sequence-family-disjoint',
    '  or new external test set, with the representation and selection policy fixed first.',
    '', '## Reproducibility', '',
    'See [the run report](abp_gn_fft1024_seed7/report.json),',
    '[all metrics](abp_gn_fft1024_seed7/metrics.csv),',
    '[predictions](abp_gn_fft1024_seed7/test_predictions.csv),',
    'and [the run instructions](README.md). Exact input files and source files are hashed;',
    'split indices and ordered feature names are saved. The benchmark checks the',
    'legacy calculation against the upstream source in its regression tests.', '',
    '![Correction audit](abp_gn_fft1024_seed7/correction_audit.png)', '']
(ROOT/'ABP_GN_RESULTS.md').write_text('\n'.join(lines))
print((ROOT/'ABP_GN_RESULTS.md').read_text())
