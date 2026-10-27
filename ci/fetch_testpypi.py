"""Download and verify the exact wheel and sdist from TestPyPI or PyPI."""

import argparse
import hashlib
import json
from pathlib import Path
import time
from urllib.error import HTTPError
from urllib.parse import urlparse
from urllib.request import urlopen


def main():
    """Check uploaded artifacts against the candidate before testing installs."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--version', required=True)
    parser.add_argument('--candidate-dir', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--index', choices=['testpypi', 'pypi'], default='testpypi')
    args = parser.parse_args()
    domain = 'test.pypi.org' if args.index == 'testpypi' else 'pypi.org'
    file_host = 'test-files.pythonhosted.org' if args.index == 'testpypi' else 'files.pythonhosted.org'
    url = f'https://{domain}/pypi/protspect/{args.version}/json'
    for attempt in range(10):
        try:
            with urlopen(url, timeout=30) as response:
                metadata = json.load(response)
            break
        except HTTPError as error:
            if error.code != 404 or attempt == 9:
                raise
            time.sleep(5)
    assert metadata['info']['name'].lower() == 'protspect'
    assert metadata['info']['version'] == args.version
    candidates = {p.name: p for p in args.candidate_dir.iterdir()
                  if p.name.endswith(('.whl', '.tar.gz'))}
    files = {item['filename']: item for item in metadata['urls']}
    if len(candidates) != 2 or files.keys() != candidates.keys():
        raise ValueError('TestPyPI files do not match the two candidate artifacts')
    args.output_dir.mkdir(parents=True, exist_ok=False)
    for name, item in files.items():
        parsed = urlparse(item['url'])
        if parsed.scheme != 'https' or parsed.hostname != file_host:
            raise ValueError('Unexpected package-index artifact host')
        expected = hashlib.sha256(candidates[name].read_bytes()).hexdigest()
        if item['digests']['sha256'] != expected:
            raise ValueError(f'Uploaded digest differs from candidate: {name}')
        with urlopen(item['url'], timeout=60) as response:
            content = response.read()
        if hashlib.sha256(content).hexdigest() != expected:
            raise ValueError(f'Download digest mismatch: {name}')
        (args.output_dir/name).write_bytes(content)
        print(f'Verified {args.index} download: {name} sha256={expected}')
    (args.output_dir/f'{args.index}-metadata.json').write_text(
        json.dumps(metadata, indent=2)+'\n', encoding='utf-8')


if __name__ == '__main__':
    main()
