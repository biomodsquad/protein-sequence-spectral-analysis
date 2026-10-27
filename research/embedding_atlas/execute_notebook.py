"""Execute the atlas notebook with the calling Python in a workspace-local kernel."""

import argparse
import json
import os
from pathlib import Path
import sys

import nbformat
from nbclient import NotebookClient


def main():
    """Keep kernel caches local and preserve earlier executed notebooks."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('Choose a new output notebook to preserve the earlier execution')
    root = Path(__file__).resolve().parents[2]
    cache = root/'.aaindex'/'jupyter'
    share = cache/'share'/'jupyter'
    kernel = share/'kernels'/'protspect-atlas'
    kernel.mkdir(parents=True,exist_ok=True)
    (kernel/'kernel.json').write_text(json.dumps({
        'argv':[sys.executable,'-m','ipykernel_launcher','-f','{connection_file}'],
        'display_name':'Protspect atlas (current Python)','language':'python',
    }),encoding='utf-8')
    environment = os.environ.copy()
    environment['JUPYTER_PATH'] = str(share)+os.pathsep+environment.get('JUPYTER_PATH','')
    for name,folder in [('MPLCONFIGDIR','matplotlib'),('NUMBA_CACHE_DIR','numba'),
                        ('IPYTHONDIR','ipython'),('JUPYTER_CONFIG_DIR','config')]:
        path = cache/folder
        path.mkdir(parents=True,exist_ok=True)
        environment[name] = str(path)
    environment.update(NUMBA_NUM_THREADS='1',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',
                       PYDEVD_DISABLE_FILE_VALIDATION='1')
    notebook = nbformat.read(Path(__file__).with_name('Embedding_atlas.ipynb'),as_version=4)
    # Kernel discovery occurs in the parent process before the subprocess inherits env.
    old = os.environ.get('JUPYTER_PATH')
    os.environ['JUPYTER_PATH'] = environment['JUPYTER_PATH']
    try:
        client = NotebookClient(notebook,timeout=900,kernel_name='protspect-atlas')
        client.execute(cwd=str(root),env=environment)
    finally:
        if old is None:
            os.environ.pop('JUPYTER_PATH',None)
        else:
            os.environ['JUPYTER_PATH'] = old
    args.output.parent.mkdir(parents=True,exist_ok=True)
    nbformat.write(notebook,args.output)
    print('Executed:',args.output)
    print('Code cells:',sum(c.cell_type == 'code' for c in notebook.cells))


if __name__ == '__main__':
    main()
