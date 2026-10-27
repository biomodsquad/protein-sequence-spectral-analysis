"""Execute the three self-contained tutorial notebooks in a clean working directory."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
from time import perf_counter

import nbformat
from nbclient import NotebookClient


def main():
    """Verify real notebook execution using this Python, without external datasets."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    notebooks = sorted((root / 'protspect' / 'examples').glob('*.ipynb'))
    if len(notebooks) != 3:
        raise ValueError('Expected the three curated offline tutorials')
    output = args.output.resolve()
    if output.exists():
        raise FileExistsError('Choose a fresh output directory to preserve earlier executions')
    output.mkdir(parents=True)
    report = []
    with TemporaryDirectory(prefix='protspect-tutorials-') as temporary:
        workspace = Path(temporary)
        share = workspace / 'share' / 'jupyter'
        kernel = share / 'kernels' / 'protspect-docs'
        kernel.mkdir(parents=True)
        (kernel / 'kernel.json').write_text(json.dumps({
            'argv': [sys.executable, '-m', 'ipykernel_launcher', '-f', '{connection_file}'],
            'display_name': 'Protspect tutorials', 'language': 'python'}))
        environment = os.environ.copy()
        environment['JUPYTER_PATH'] = str(share) + os.pathsep + environment.get('JUPYTER_PATH', '')
        environment['PYTHONPATH'] = str(root) + os.pathsep + environment.get('PYTHONPATH', '')
        environment.update(PYDEVD_DISABLE_FILE_VALIDATION='1', OMP_NUM_THREADS='1',
                           OPENBLAS_NUM_THREADS='1', NUMBA_NUM_THREADS='1')
        for variable, folder in [('MPLCONFIGDIR', 'matplotlib'), ('IPYTHONDIR', 'ipython'),
                                 ('JUPYTER_CONFIG_DIR', 'config')]:
            directory = workspace / folder
            directory.mkdir()
            environment[variable] = str(directory)
        old = os.environ.get('JUPYTER_PATH')
        os.environ['JUPYTER_PATH'] = environment['JUPYTER_PATH']
        try:
            for path in notebooks:
                notebook = nbformat.read(path, as_version=4)
                nbformat.validate(notebook)
                if any(c.get('outputs') for c in notebook.cells):
                    raise ValueError(f'Source tutorial must have clean outputs: {path.name}')
                start = perf_counter()
                NotebookClient(notebook, timeout=600, kernel_name='protspect-docs').execute(
                    cwd=str(workspace), env=environment)
                cells = [c for c in notebook.cells if c.cell_type == 'code']
                if any(c.execution_count is None for c in cells):
                    raise RuntimeError('A code cell did not execute')
                if any(o.output_type == 'error' for c in cells for o in c.outputs):
                    raise RuntimeError('An error output was recorded')
                nbformat.write(notebook, output / path.name)
                report.append({'notebook': path.name, 'code_cells': len(cells),
                    'seconds': round(perf_counter()-start, 3),
                    'source_sha256': hashlib.sha256(path.read_bytes()).hexdigest()})
                print('Executed', path.name, flush=True)
        finally:
            if old is None:
                os.environ.pop('JUPYTER_PATH', None)
            else:
                os.environ['JUPYTER_PATH'] = old
    (output / 'report.json').write_text(json.dumps(report, indent=2)+'\n')
    print('All tutorials passed; results:', output)


if __name__ == '__main__':
    main()
