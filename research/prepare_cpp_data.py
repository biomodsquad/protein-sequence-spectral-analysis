"""Prepare complete, hash-verified MLCPP 2.0 Layer1 downloads for the example."""

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
from urllib.request import urlopen

import pandas as pd

MANIFEST = Path(__file__).resolve().with_name('mlcpp2-source.json')


def parse_fasta(text):
    """Preserve official FASTA identifiers, order, sequences, and labels."""
    rows, identifier, chunks = [], None, []
    for line in text.splitlines() + ['>END']:
        if line.startswith('>'):
            if identifier is not None:
                label = identifier.split('_')[0]
                sequence = ''.join(chunks)
                if label not in ('Positive', 'Negative') or not sequence or set(sequence) - set('ACDEFGHIKLMNPQRSTVWY'):
                    raise ValueError('invalid official FASTA record')
                rows.append({'SourceID': identifier, 'Sequence': sequence, 'Class': label})
            identifier, chunks = line[1:].strip(), []
        elif line.strip():
            if identifier is None:
                raise ValueError('sequence before FASTA identifier')
            chunks.append(line.strip())
    if not rows:
        raise ValueError('empty FASTA data')
    table = pd.DataFrame(rows)
    if table.SourceID.duplicated().any() or table.Sequence.duplicated().any():
        raise ValueError('duplicate official records require explicit resolution')
    return table


def prepare(output_dir, source_dir=None):
    """Download or read cached official files and verify bytes before writing data.

    The committed manifest pins the reviewed website snapshot. Changed downloads
    fail the hash check, requiring an explicit provenance review. No records are
    filtered or reordered. Source files and generated CSVs stay in a local cache.
    """
    output = Path(output_dir)
    if output.exists():
        raise FileExistsError('choose a fresh data directory to preserve the existing snapshot')
    manifest = json.loads(MANIFEST.read_text())
    prepared, provenance = {}, {'source_page': manifest['source_page'],
                               'snapshot_date': manifest['snapshot_date'], 'datasets': {}}
    for name, spec in manifest['datasets'].items():
        if source_dir is None:
            with urlopen(spec['url'], timeout=30) as response:
                raw = response.read()
        else:
            raw = (Path(source_dir) / f'mlcpp2_{name}.fasta').read_bytes()
        digest = hashlib.sha256(raw).hexdigest()
        if digest != spec['sha256']:
            raise ValueError(f'{name} SHA-256 differs from the reviewed official snapshot')
        table = parse_fasta(raw.decode('utf-8-sig'))
        if len(table) != spec['rows'] or dict(Counter(table.Class)) != spec['classes']:
            raise ValueError(f'{name} record counts differ from the official manifest')
        csv = table.to_csv(index=False, lineterminator='\n').encode()
        prepared[name] = (raw, csv)
        provenance['datasets'][name] = {**spec, 'csv_file': f'{name}.csv',
                                      'csv_sha256': hashlib.sha256(csv).hexdigest()}
    output.mkdir(parents=True)
    for name, (raw, csv) in prepared.items():
        (output / f'mlcpp2_{name}.fasta').write_bytes(raw)
        (output / f'{name}.csv').write_bytes(csv)
    (output/'source.json').write_text(json.dumps(provenance, indent=2) + '\n')
    return provenance


def main():
    """Prepare the reviewed website snapshot, optionally using verified cached bytes."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=Path('research/data/mlcpp2'))
    parser.add_argument('--source-dir', type=Path, help='cache with mlcpp2_training/independent.fasta')
    args = parser.parse_args()
    print(json.dumps(prepare(args.output_dir, args.source_dir), indent=2))


if __name__ == '__main__':
    main()
