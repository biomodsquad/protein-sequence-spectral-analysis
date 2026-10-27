"""Execute the AMP transform notebook using the calling Python environment."""

import argparse
import json
import os
from pathlib import Path
import sys

import nbformat
from nbclient import NotebookClient


def main():
    """Run the full benchmark and preserve an executed copy with its outputs."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('Choose a new output path to preserve earlier results')
    root = Path(__file__).resolve().parents[1]
    cache = root / '.aaindex' / 'amp_notebook'
    share = cache / 'share' / 'jupyter'
    kernel = share / 'kernels' / 'protspect-amp'
    kernel.mkdir(parents=True, exist_ok=True)
    (kernel / 'kernel.json').write_text(json.dumps({
        'argv': [sys.executable, '-m', 'ipykernel_launcher', '-f', '{connection_file}'],
        'display_name': 'Protspect AMP (current Python)', 'language': 'python'}))
    environment = os.environ.copy()
    environment['JUPYTER_PATH'] = str(share) + os.pathsep + environment.get('JUPYTER_PATH', '')
    for variable, folder in [('MPLCONFIGDIR', 'matplotlib'), ('IPYTHONDIR', 'ipython'),
                             ('JUPYTER_CONFIG_DIR', 'config')]:
        path = cache / folder
        path.mkdir(parents=True, exist_ok=True)
        environment[variable] = str(path)
    environment['PYDEVD_DISABLE_FILE_VALIDATION'] = '1'
    notebook = nbformat.read(root / 'research' / 'AMP_transform_MLP.ipynb', as_version=4)
    old = os.environ.get('JUPYTER_PATH')
    os.environ['JUPYTER_PATH'] = environment['JUPYTER_PATH']
    try:
        def progress(cell, cell_index, execute_reply):
            print(f'Completed cell {cell_index + 1}/{len(notebook.cells)}', flush=True)
            for output in cell.get('outputs', []):
                if output.get('output_type') == 'stream':
                    print(output.get('text', ''), end='', flush=True)
        NotebookClient(notebook, timeout=3600, kernel_name='protspect-amp',
                       on_cell_executed=progress).execute(cwd=str(root), env=environment)
    finally:
        if old is None:
            os.environ.pop('JUPYTER_PATH', None)
        else:
            os.environ['JUPYTER_PATH'] = old
    args.output.parent.mkdir(parents=True, exist_ok=True)
    nbformat.write(notebook, args.output)
    print('Executed notebook:', args.output)


if __name__ == '__main__':
    main()
