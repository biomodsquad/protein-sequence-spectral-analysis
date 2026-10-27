"""Generate matrix-profile features from an explicitly supplied local AAindex file."""

import argparse
import json
from pathlib import Path

from protspect import SequenceSet
from protspect.aaindex import AAIndex2, AAIndex3


def main():
    """Require the matrix section and profile axis; optionally export provenance."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('database', type=Path)
    parser.add_argument('accessions', nargs='+')
    parser.add_argument('--section', type=int, choices=(2, 3), required=True)
    parser.add_argument('--axis', choices=('rows', 'columns'), required=True)
    parser.add_argument('--expected-sha256')
    parser.add_argument('--metadata-output', type=Path)
    args = parser.parse_args()
    catalogue = AAIndex2 if args.section == 2 else AAIndex3
    db = catalogue.from_file(args.database, expected_sha256=args.expected_sha256)
    features = SequenceSet(embedding=db.embedding(args.accessions, axis=args.axis),
                           compute_cross_spectra=False)
    features.from_sequence_list(['ACDE', 'ACDF'], spectrum_length=8)
    print(features.export_datatable())
    if args.metadata_output:
        args.metadata_output.write_text(json.dumps(db.metadata(args.accessions, axis=args.axis),
                                                  indent=2, allow_nan=False), encoding='utf-8')


if __name__ == '__main__':
    main()
