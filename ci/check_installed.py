"""Test wheel and sdist installs in fresh venvs outside the source tree, on any OS."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import venv

ROOT = Path(__file__).resolve().parents[1]

CORE_SMOKE = '''
import importlib.util
import sys
from pathlib import Path
import protspect
from protspect import ComplexSet, SequenceFFT, SequenceSet
from FFTpro.features import SequenceFFT as LegacyFFT
assert LegacyFFT is SequenceFFT
assert 'sklearn' not in sys.modules
assert importlib.util.find_spec('sklearn') is None
assert Path(protspect.__file__).resolve().is_relative_to(Path(sys.prefix).resolve())
s = SequenceSet(normalization='length')
s.from_sequence_list(['ACDE'], spectrum_length=8)
assert s.export_datatable().shape == (1, 80)
p = ComplexSet(normalization='energy')
p.from_sequence_list(['ACDE'], ['ACDF'], spectrum_length=8)
assert p.export_datatable().shape == (1, 1600)
try:
    import protspect.sklearn
except ImportError as error:
    assert 'Protspect[ml]' in str(error)
else:
    raise AssertionError('optional sklearn dependency unexpectedly available')
print('Core install, legacy API, normalization, paired chains, and optional imports passed.')
'''

ENVIRONMENT = '''
import importlib.metadata as metadata
import json, platform, sys
print(json.dumps({'python': platform.python_version(), 'implementation': platform.python_implementation(),
    'platform': platform.system(), 'machine': platform.machine(),
    'packages': {d.metadata['Name']: d.version for d in metadata.distributions()}}, sort_keys=True))
'''

PLOT_SMOKE = '''
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
fig, ax = plt.subplots()
ax.plot([0, 1], [0, 1])
fig.savefig('optional-validation-smoke.png')
plt.close(fig)
'''


def command(args, cwd, env, log):
    """Run a check with an exact argument vector and retain useful failure logs."""
    result = subprocess.run([str(arg) for arg in args], cwd=cwd, env=env,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    with log.open('a', encoding='utf-8') as stream:
        stream.write('\nCOMMAND ' + repr([str(arg) for arg in args]) + '\n' + result.stdout)
    print(result.stdout, end='', flush=True)
    if result.returncode:
        raise subprocess.CalledProcessError(result.returncode, args)
    return result.stdout


def check_artifact(artifact, output, constraints=None):
    """Install core, run its tests, then add validation extras and run the full suite."""
    artifact, output = Path(artifact).resolve(), Path(output).resolve()
    output.mkdir(parents=True, exist_ok=False)
    env = os.environ.copy()
    for name in ('PYTHONPATH', 'PYTHONHOME', 'VIRTUAL_ENV', 'PIP_TARGET', 'PIP_PREFIX'):
        env.pop(name, None)
    env.update(PYTHONNOUSERSITE='1', PIP_DISABLE_PIP_VERSION_CHECK='1',
               MPLCONFIGDIR=str(output/'matplotlib'), PYTHONUTF8='1')
    if constraints:
        env['PIP_CONSTRAINT'] = Path(constraints).resolve().as_uri()
    log = output/'checks.log'
    with tempfile.TemporaryDirectory(prefix='protspect-installed-') as temporary:
        work = Path(temporary)
        environment = work/'venv'
        venv.EnvBuilder(with_pip=True, system_site_packages=False).create(environment)
        python = environment/('Scripts/python.exe' if os.name == 'nt' else 'bin/python')
        # The copied tests can see their example helpers but cannot import package
        # source. No .git exists here, so the checkout-only historical-source test skips.
        suite = work/'suite'
        shutil.copytree(ROOT/'tests', suite/'tests', ignore=shutil.ignore_patterns('__pycache__'))
        shutil.copytree(ROOT/'examples', suite/'examples', ignore=shutil.ignore_patterns('__pycache__'))
        command([python, '-m', 'pip', 'install', artifact, 'pytest>=7'], work, env, log)
        command([python, '-m', 'pip', 'check'], work, env, log)
        core_env = json.loads(command([python, '-c', ENVIRONMENT], work, env, log))
        command([python, '-c', CORE_SMOKE], work, env, log)
        command([python, '-m', 'pytest', '-q', '-ra', suite/'tests'], work, env, log)
        command([python, ROOT/'examples/basic_features.py'], work, env, log)
        command([python, '-m', 'pip', 'install', 'Protspect[validation] @ '+artifact.as_uri()], work, env, log)
        command([python, '-m', 'pip', 'check'], work, env, log)
        optional_env = json.loads(command([python, '-c', ENVIRONMENT], work, env, log))
        command([python, '-m', 'pytest', '-q', '-ra', suite/'tests'], work, env, log)
        command([python, ROOT/'examples/normalized_pipeline.py'], work, env, log)
        command([python, '-c', PLOT_SMOKE], work, env, log)
    result = {'artifact': artifact.name, 'sha256': hashlib.sha256(artifact.read_bytes()).hexdigest(),
              'constraints': None if constraints is None else Path(constraints).name,
              'core': core_env, 'validation': optional_env, 'checks_passed': True,
              'test_scope': 'installed package; checkout-only git comparison skipped; optional modules skipped in core-only stage'}
    (output/'environment.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    return result


def main():
    """Select an interpreter by running this helper with that Python executable."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dist-dir', type=Path, default=Path('dist'))
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--constraints', type=Path)
    args = parser.parse_args()
    artifacts = sorted(args.dist_dir.glob('*.whl')) + sorted(args.dist_dir.glob('*.tar.gz'))
    if len(artifacts) != 2:
        raise ValueError('dist-dir must contain exactly one wheel and one source archive')
    for artifact in artifacts:
        label = 'wheel' if artifact.suffix == '.whl' else 'sdist'
        check_artifact(artifact, args.output_dir/label, args.constraints)
    print('All fresh installed-distribution checks passed on', sys.version.split()[0])


if __name__ == '__main__':
    main()
