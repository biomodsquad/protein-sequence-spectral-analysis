"""Export the current executed notebook's seven Plotly figures for offline docs.

Run execute_notebook --plotly-json first. Figures are copied from its outputs;
this exporter does not recompute layouts or require AAindex downloads.
"""

import argparse
import hashlib
import json
from pathlib import Path

import nbformat
import plotly
import plotly.graph_objects as go
import plotly.io as pio
from plotly.offline import get_plotlyjs

FIGURES = [
    ('vector-pca', 'Individual-channel PCA'),
    ('vector-umap', 'Individual channels: UMAP'),
    ('whole-pca', 'Whole residue-geometry PCA'),
    ('whole-umap', 'Whole residue geometry: UMAP'),
    ('cka-heatmap', 'Whole representations: original-space CKA with dendrograms'),
    ('channel-matrix', 'BLOSUM100: 20/20 channels displayed'),
    ('matrix-pca', 'Signed aligned matrix-score PCA'),
]


def export(notebook_path, output):
    """Validate expected figures and save a shared JS bundle plus provenance."""
    notebook = nbformat.read(notebook_path, as_version=4)
    figures = [item['data']['application/vnd.plotly.v1+json']
               for cell in notebook.cells if cell.cell_type == 'code'
               for item in cell.get('outputs', [])
               if 'application/vnd.plotly.v1+json' in item.get('data', {})]
    if len(figures) != len(FIGURES):
        raise ValueError(f'Expected {len(FIGURES)} Plotly figures, found {len(figures)}')
    for (_, title), figure in zip(FIGURES, figures):
        if figure.get('layout', {}).get('title', {}).get('text') != title:
            raise ValueError(f'Unexpected notebook figure; expected {title!r}')
    output.mkdir(parents=True, exist_ok=True)
    (output/'plotly.min.js').write_text(get_plotlyjs(), encoding='utf-8')
    records = []
    for (slug, title), figure in zip(FIGURES, figures):
        path = output/f'{slug}.html'
        pio.write_html(go.Figure(figure), path, include_plotlyjs='plotly.min.js',
                       full_html=True, auto_open=False, div_id=slug,
                       config={'responsive': True, 'displaylogo': False})
        records.append({'file': path.name, 'title': title,
                        'sha256': hashlib.sha256(path.read_bytes()).hexdigest()})
    root = Path(__file__).resolve().parents[2]
    sources = ['research/embedding_atlas/Embedding_atlas.ipynb',
               'research/embedding_atlas/analyze.py', 'research/embedding_atlas/explore.py',
               'research/embedding_atlas/execute_notebook.py',
               'research/embedding_atlas/export_docs.py']
    report = {'plotly_version': plotly.__version__, 'figures': records,
              'executed_notebook_sha256': hashlib.sha256(notebook_path.read_bytes()).hexdigest(),
              'sources': {name: hashlib.sha256((root/name).read_bytes()).hexdigest()
                          for name in sources},
              'aaindex_sha256': {f'aaindex{i}': hashlib.sha256(
                  (root/'.aaindex'/f'aaindex{i}').read_bytes()).hexdigest() for i in (1, 2, 3)}}
    (output/'manifest.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print(f'Exported {len(figures)} interactive figures to {output}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('notebook', type=Path)
    parser.add_argument('--output', type=Path, default=Path('docs/_static/embedding-atlas'))
    args = parser.parse_args()
    export(args.notebook, args.output)
