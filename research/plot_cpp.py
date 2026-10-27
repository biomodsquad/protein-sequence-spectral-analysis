"""Plot recorded decisions and feature attributions; never fit or select a model."""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from sklearn.metrics import PrecisionRecallDisplay, RocCurveDisplay


def plot_results(directory):
    """Save ROC/PR and optional MISTIC attribution figures from run artifacts."""
    plt.switch_backend('Agg')
    directory = Path(directory)
    data = pd.read_csv(directory/'test_predictions.csv')
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    for name, label in [('rbf_svc', 'RBF SVC'), ('mistic', 'MISTIC compact search')]:
        if name + '_decision' not in data:
            continue
        RocCurveDisplay.from_predictions(data.Class, data[name + '_decision'], name=label, ax=axes[0])
        PrecisionRecallDisplay.from_predictions(data.Class, data[name + '_decision'], name=label, ax=axes[1])
    axes[0].plot([0, 1], [0, 1], '--', color='gray', linewidth=1)
    axes[1].axhline(data.Class.mean(), linestyle='--', color='gray', linewidth=1,
                   label=f'Prevalence {data.Class.mean():.1%}')
    axes[1].legend(fontsize=8)
    fig.suptitle('CPP example · exact train/test duplicates excluded')
    fig.savefig(directory/'evaluation.png', dpi=170)
    plt.close(fig)
    if (directory/'attribution_summary.csv').exists():
        importance = pd.read_csv(directory/'attribution_summary.csv').head(12).iloc[::-1]
        fig, ax = plt.subplots(figsize=(8, 3.6), constrained_layout=True)
        ax.barh(importance.feature, importance.mean_absolute_integrated_gradient, color='#337e92')
        ax.set_xlabel('Mean absolute integrated gradient · 32 fixed example rows')
        ax.set_title('MISTIC spectral-feature explanation · zero standardized-feature baseline')
        ax.spines[['top', 'right']].set_visible(False)
        fig.savefig(directory/'attributions.png', dpi=170)
        plt.close(fig)
