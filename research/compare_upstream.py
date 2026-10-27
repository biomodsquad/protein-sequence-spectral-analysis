"""Compare local CPP CSVs with labeled FASTA downloads without redistributing sequences."""

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

import pandas as pd


def fasta_records(path):
    """Parse wrapped FASTA records whose identifiers start Positive_ or Negative_."""
    result, identifier, sequence = [], None, []
    for line in Path(path).read_text().splitlines():
        if line.startswith('>'):
            if identifier is not None:
                result.append((''.join(sequence), identifier.split('_')[0]))
            identifier, sequence = line[1:].strip(), []
        elif line.strip():
            if identifier is None:
                raise ValueError('sequence before FASTA header')
            sequence.append(line.strip())
    if identifier is not None:
        result.append((''.join(sequence), identifier.split('_')[0]))
    if not result or any(not sequence or label not in ('Positive', 'Negative') for sequence, label in result):
        raise ValueError('expected nonempty Positive_/Negative_ FASTA records')
    return result


def compare(csv, fasta, url):
    """Compare labeled multisets, preserving duplicate multiplicity and exact source hashes."""
    records = fasta_records(fasta)
    frame = pd.read_csv(csv)
    official, repository = Counter(records), Counter(zip(frame.Sequence, frame.Class))
    omitted_classes = Counter()
    for (_, label), count in (official - repository).items():
        omitted_classes[label] += count
    def digest(path):
        return hashlib.sha256(Path(path).read_bytes()).hexdigest()
    return {'official_url': url, 'official_sha256': digest(fasta), 'repository_sha256': digest(csv),
            'official_rows': len(records), 'official_classes': dict(Counter(label for _, label in records)),
            'repository_rows': len(frame), 'matching_labeled_records': sum((repository & official).values()),
            'repository_only': sum((repository - official).values()),
            'official_only': sum((official - repository).values()),
            'official_only_classes': dict(omitted_classes), 'same_labeled_multiset': repository == official}


def main():
    """Audit the two official downloads against the repository CSVs."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--training-fasta', type=Path, required=True)
    parser.add_argument('--test-fasta', type=Path, required=True)
    parser.add_argument('--csv-dir', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    base = 'https://balalab-skku.org/mlcpp2/download/'
    report = {'source_page': base, 'checked_date': '2026-10-27', 'comparisons': {}}
    for filename, fasta, remote in [('CPP_training.csv', args.training_fasta, 'Layer1_training.txt'),
                                    ('CPP_independent.csv', args.test_fasta, 'Layer1_Independent.txt')]:
        report['comparisons'][filename] = compare(args.csv_dir / filename, fasta, base + remote + '/txt/')
    args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
