"""Generate reproducible documentation figures from actual Protspect outputs."""

from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import numpy as np

from protspect import transform_sequence
from protspect.embeddings import b100, kidera_factors

OUTPUT = Path(__file__).parent / '_static' / 'figures'


def save(figure, name):
    """Save scalable site figures and raster copies for visual review."""
    OUTPUT.mkdir(parents=True, exist_ok=True)
    figure.savefig(OUTPUT / f'{name}.svg', bbox_inches='tight',
                   metadata={'Date': None, 'Creator': 'Protspect tutorial figure generator'})
    figure.savefig(OUTPUT / f'{name}.png', dpi=150, bbox_inches='tight')
    plt.close(figure)


def workflow():
    """Show the shared processing stages, independently of the provided diagram."""
    figure, axis = plt.subplots(figsize=(12, 2.2))
    axis.set(xlim=(0, 12), ylim=(0, 2))
    axis.axis('off')
    labels = ['Sequence\nN residues', 'Embedding\nN × C signals',
              'Transform\nalong position', 'Features\nselect coefficients', 'Model\nfit on training']
    colors = ['#e9f2fb', '#e4f5f2', '#edf2fd', '#fff3dc', '#e9f2fb']
    for index, (label, color) in enumerate(zip(labels, colors)):
        x = .1 + index * 2.4
        axis.add_patch(FancyBboxPatch((x, .55), 2.05, .95,
            boxstyle='round,pad=.08', facecolor=color, edgecolor='#7391a9'))
        axis.text(x + 1.025, 1.025, label, ha='center', va='center', fontsize=12, color='#173d5c')
        if index < 4:
            axis.annotate('', xy=(x+2.34, 1), xytext=(x+2.15, 1),
                          arrowprops={'arrowstyle': '->', 'color': '#567487'})
    axis.text(6, .1, 'Preprocess real residues before zero-padding · Keep learned statistics inside each training fold',
              ha='center', fontsize=11, color='#34566b')
    save(figure, 'workflow')


def representations():
    """Plot raw residue/channel score matrices with separate scales."""
    alphabet = 'ACDEFGHIKLMNPQRSTVWY'
    figure, axes = plt.subplots(1, 2, figsize=(11, 5.3), constrained_layout=True)
    for axis, (name, embedding) in zip(axes, [('BLOSUM100', b100), ('Kidera', kidera_factors())]):
        channels = list(embedding)
        values = np.array([[embedding[c][aa] for c in channels] for aa in alphabet])
        peak = abs(values).max()
        image = axis.imshow(values, cmap='RdBu_r', vmin=-peak, vmax=peak, aspect='auto')
        axis.set_xticks(np.arange(len(channels)), channels, rotation=90, fontsize=8)
        axis.set_yticks(np.arange(20), list(alphabet))
        axis.set(xlabel='Embedding channel', ylabel='Residue', title=f'{name}: 20 × {len(channels)}')
        figure.colorbar(image, ax=axis, label='Raw channel score', shrink=.85)
    save(figure, 'representations')


def transforms():
    """Compare real coefficients while keeping each transform's axis labels."""
    options = {'fourier': {}, 'z': {'radius': 1.05}, 'mellin': {'mellin_c': 0.},
               'wavelet': {'wavelet': 'haar', 'level': 3, 'mode': 'periodization'}}
    figure, axes = plt.subplots(2, 2, figsize=(11, 6.4), constrained_layout=True)
    for axis, (method, settings) in zip(axes.flat, options.items()):
        result = transform_sequence('IKLFYRACDEFGHIK', kidera_factors(), method=method,
                                    length=64, normalization='energy', **settings)
        values = result.coefficients['KF1'].real
        if method == 'wavelet':
            x, label = np.arange(len(values)), 'Concatenated band coefficient'
        elif method == 'mellin':
            x, label = result.coordinates['log_frequency'], 'Log-frequency ν'
        else:
            x, label = result.coordinates['frequency'], 'Cycles per residue'
        axis.plot(x, values, color='#205b87', linewidth=1.5)
        axis.axhline(0, color='#8999a5', linewidth=.6)
        axis.set(title=method.capitalize(), xlabel=label, ylabel='Real coefficient')
        axis.grid(alpha=.15)
    save(figure, 'transforms')


def main():
    """Rebuild all generated figures; user-supplied artwork is never modified."""
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'svg.fonttype': 'none'})
    workflow()
    representations()
    transforms()
    print('Generated figures in', OUTPUT)


if __name__ == '__main__':
    main()
