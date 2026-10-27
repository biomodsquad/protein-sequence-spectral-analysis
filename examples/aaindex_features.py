"""Generate features from an explicitly provided local AAindex1 database."""

import argparse
import json
from pathlib import Path

from protspect import SequenceSet
from protspect.aaindex import AAIndex


def main():
    """Load offline, select property channels, and optionally save source provenance."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('database', type=Path)
    parser.add_argument('--accessions', nargs='+', default=['KYTJ820101', 'FAUJ880103'])
    parser.add_argument('--expected-sha256')
    parser.add_argument('--metadata-output', type=Path)
    args = parser.parse_args()
    database = AAIndex.from_file(args.database, expected_sha256=args.expected_sha256)
    features = SequenceSet(embedding=database.embedding(args.accessions),
                           compute_cross_spectra=False)
    features.from_sequence_list(['ACDE', 'ACDF'], spectrum_length=8)
    print(features.export_datatable())
    if args.metadata_output:
        args.metadata_output.write_text(json.dumps(database.metadata(args.accessions), indent=2),
                                        encoding='utf-8')


if __name__ == '__main__':
    main()
