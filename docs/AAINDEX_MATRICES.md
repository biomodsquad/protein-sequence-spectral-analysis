# AAindex2 and AAindex3 matrix profiles

The development version supports all three AAindex file formats. These
additions are not in published PyPI 0.1.0a1. BLOSUM100 remains the default;
AAindex matrix conversion does not change Fourier calculations.

## Explicit downloads and offline catalogues

```python
from protspect.aaindex import (
    AAIndex2, AAIndex3, download_aaindex2, download_aaindex3,
)

# Run once; existing files are never overwritten.
download_aaindex2(".aaindex/aaindex2")
download_aaindex3(".aaindex/aaindex3")

# Reuse locally without network access.
substitutions = AAIndex2.from_file(".aaindex/aaindex2")
contacts = AAIndex3.from_file(".aaindex/aaindex3")
```

All functions accept `expected_sha256` for raw-file pinning. Download helpers
validate records and checksum before writing. No database is bundled. The two
matrix sections share a file format: choose the correct catalogue explicitly;
the class labels provenance and does not infer biological meaning from numbers.
For scalar indices, keep using `AAIndex` as described in [AAindex1](AAINDEX.md).

Search, iteration, lookup, version/source labels, and checksum behavior match
the scalar catalogue. `search(..., complete_only=True)` tests the entire native
matrix. Individual entries expose `rows`, `columns`, `shape`, `storage`,
`values[row][column]`, `missing_pairs`, `complete`, and `canonical_complete`.
Entries and nested value mappings are read-only.

## Choose profile direction explicitly

```python
from protspect import SequenceSet

embedding = substitutions.embedding("HENS920102", axis="rows")  # BLOSUM62
features = SequenceSet(embedding=embedding, compute_cross_spectra=False)
features.from_sequence_list(["ACDE", "ACDF"], spectrum_length=8)
table = features.export_datatable()
```

A selected matrix produces **20 channels**, each covering the 20 canonical
residues in `ARNDCQEGHILKMFPSTWYV` order. Both orientations are available:

| Selection | Channel | Value assigned to residue R |
| --- | --- | --- |
| `axis="rows"` | `<accession>_rows_A` | `M[A, R]` |
| `axis="columns"` | `<accession>_columns_A` | `M[R, A]` |

The `axis` keyword is required. It affects asymmetric matrices; symmetric
matrices give equivalent numbers with different orientation-labelled names.
Channel names include accession and axis so multiple matrices remain distinct.
Requested accession order is preserved. Returned ordinary dictionaries work
with `SequenceSet`, `ComplexSet`, and `SpectralTransformer`, including cloning
and pickle serialization. Raw scores are retained without normalization.

The adapter explicitly selects the canonical 20-by-20 block. Gaps, J/O, and
other native symbols remain inspectable in the original matrix and provenance,
but are excluded from these profile channels. Missing canonical residues or
cells raise `ValueError`; no imputation, averaging, or alphabet expansion is
performed. Missing values outside the canonical block do not invalidate that
block. Unknown sequence symbols retain the normal core validation behavior.

## Storage and direction

The parser reads the M header's separate row and column alphabets. Full
matrices preserve every row/column value, even if rectangular or asymmetric.
Only explicitly lower-triangular storage with identical ordered axes is
reflected across the diagonal. Upper-triangular or malformed layouts are
rejected rather than guessed. Both `NA` and a standalone `-` score are retained
as missing; `-` in an axis is instead a gap symbol. Underscore accessions such
as `QU_C930101` are supported.

Preserve the source header and comments when interpreting values. Protspect
uses the declared numeric axes and does not silently reinterpret or repair
possible inconsistencies between source comments and those axes.

## AAindex3 interpretation

```python
embedding = contacts.embedding("TANS760101", axis="columns")
```

Contact-potential profiles assign twenty reference-residue scores to each
sequence residue. Their Fourier spectra are **sequence-derived features**;
they do not predict contacts, calculate a folded protein's contact energy, or
supply spatial information absent from the sequence. Preserve the original
source publication's sign conventions and units. No predictive advantage is
claimed without an appropriate benchmark.

## Save source and conversion provenance

```python
import json
from pathlib import Path

metadata = contacts.metadata("TANS760101", axis="columns")
Path(".aaindex/contact-profile-source.json").write_text(
    json.dumps(metadata, indent=2, allow_nan=False), encoding="utf-8"
)
```

Supply the same axis as used for embedding conversion. Metadata records section,
raw checksum, local path/source URL/version label, full native axes and values,
input storage form, original citation/comments, selected channel names, profile
orientation, and canonical-submatrix policy. Save it alongside the estimator's
`get_feature_metadata()` output. Calling `metadata(ids)` without an axis
exports source records, including incomplete ones, without claiming a derived
profile schema.

## Full-file validation

Official snapshots inspected on 2026-10-27:

| Section | Records | Lower triangle / full | Incomplete native matrices | SHA256 |
| --- | --- | --- | --- | --- |
| AAindex2 | 94 | 67 / 27 | 2 | `d3eb30c075452aa737b111a3d127358d0b669d0e847fca967b253a05c7796aa9` |
| AAindex3 | 47 | 44 / 3 | 4 | `8c9b9b0e7b7937b4674585b25f3f9370499f219b2583b1204fc9a0a5f06b1671` |

The incomplete AAindex2 entries are `MEHP950101` and `MEHP950103`; the incomplete
AAindex3 entries are `ROBB790102`, `PARB960101`, `PARB960102`, and `GODA950101`.
Both orientations were checked for every complete canonical matrix. Synthetic
fixtures exercise direction, rectangular shapes, missing cells, format errors,
metadata, downloads, and numerical spectral/estimator integration without
redistributing source databases. See [AAindex attribution](AAINDEX.md#database-attribution).

All 400 BLOSUM50 (`HENS920104`) scores match the independent
[NCBI BLOSUM50 matrix](https://ftp.ncbi.nlm.nih.gov/blast/matrices/BLOSUM50).

Source units matter: AAindex's `HENS920102` BLOSUM62 uses **1/3-bit units**,
whereas the standard NCBI BLOSUM62 file uses 1/2-bit units. Their numerical
values differ legitimately; Protspect preserves the selected source's units
and records its comments. There is no implicit conversion to another BLOSUM
variant.

See the [official matrix format](https://www.genome.jp/aaindex/aaindex_help.html).
