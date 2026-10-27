"""Execute the research notebook in its pinned kernel and retain the executed artifact."""

import argparse
import os
from pathlib import Path

import nbformat
from nbclient import NotebookClient


def main():
    """Execute from the repository root; errors fail the command and remain visible."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    if args.output_dir.exists():
        raise FileExistsError('choose a new output directory to preserve previous run artifacts')
    notebook = nbformat.read(root/'research/CPP_reproducible.ipynb', as_version=4)
    environment = os.environ.copy()
    environment['PROTSPECT_CPP_OUTPUT'] = str(args.output_dir.resolve())
    environment['MPLCONFIGDIR'] = str(root/'research/.jupyter/matplotlib')
    environment['IPYTHONDIR'] = str(root/'research/.jupyter/ipython')
    environment['JUPYTER_CONFIG_DIR'] = str(root/'research/.jupyter/config')
    environment['PYDEVD_DISABLE_FILE_VALIDATION'] = '1'
    for name in ('MPLCONFIGDIR', 'IPYTHONDIR', 'JUPYTER_CONFIG_DIR'):
        Path(environment[name]).mkdir(parents=True, exist_ok=True)
    client = NotebookClient(notebook, timeout=900, kernel_name='protspect-research')
    client.execute(cwd=str(root), env=environment)
    nbformat.write(notebook, args.output_dir/'CPP_reproducible.executed.ipynb')
    print('Executed notebook:', args.output_dir/'CPP_reproducible.executed.ipynb')


if __name__ == '__main__':
    main()
