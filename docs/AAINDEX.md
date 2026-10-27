# AAindex1 representations

Protspect 0.1.0a2 supports local AAindex1 files. This page describes scalar properties; see
[AAindex2/3 matrix profiles](AAINDEX_MATRICES.md) for substitution matrices and
contact potentials. BLOSUM100 remains the default representation.

## Download once, then work offline

No network access happens during import, local loading, fitting, or feature
extraction. Download the official database explicitly into a user-managed cache:

```python
from protspect.aaindex import AAIndex, download_aaindex1

# Run once. Existing files are never overwritten.
database = download_aaindex1(".aaindex/aaindex1")
print(database.sha256)

# Subsequent runs are completely offline. Pin the printed digest when desired.
database = AAIndex.from_file(".aaindex/aaindex1")
```

Both functions accept `expected_sha256="..."` to require an exact snapshot.
The download helper validates checksum and record structure before saving and
refuses existing destinations. Keep different snapshots at different paths.
An interrupted or invalid download does not replace an existing cache file.

AAindex1 does not declare its database release in individual records. Optional
`version` and `source_url` arguments to `from_file` are caller-supplied labels;
the raw SHA256 is always computed, never inferred from a release label. The
current official file inspected on 2026-10-27 has 566 entries, 13 with missing
values, and SHA256:

```
560f102dbee8ca9837186c4b2b93e253ea6efc58c88617234cef4e54bd61f676
```

## Find and select properties

```python
for entry in database.search("hydrophobicity", complete_only=True):
    print(entry.accession, entry.description)

entry = database["KYTJ820101"]
print(entry.values["A"])  # 1.8
print(entry.title, entry.journal)
print(entry.complete, entry.missing_residues)
```

Search matches a case-insensitive substring of the accession, description,
authors, title, journal, reference identifier, or comments. Catalogue iteration
preserves file order. Record values are read-only; missing source values remain
`None` and can be inspected without creating an embedding.

```python
from protspect import SequenceSet

ids = ["KYTJ820101", "FAUJ880103"]
embedding = database.embedding(ids)
features = SequenceSet(embedding=embedding, compute_cross_spectra=False)
features.from_sequence_list(["ACDE", "ACDF"], spectrum_length=8)
print(features.export_datatable())
```

Each accession produces one residue-property channel. Caller selection order is
preserved. Unknown IDs raise `KeyError`; duplicates or empty selections raise
`ValueError`. The parser decodes the residue-order header rather than assuming
alphabetical order. It rejects malformed records, duplicate accessions,
nonfinite values, and matrix records.

`embedding()` returns independent ordinary dictionaries, suitable for
`SequenceSet`, `ComplexSet`, and the optional `SpectralTransformer`, including
estimator cloning and serialization. The default policy rejects selected
entries containing `NA`, reporting accessions and missing residues. There is no
implicit imputation, residue filtering, averaging, or scaling. Unknown sequence
residues retain the existing Protspect validation behavior.

## Save provenance with the model

```python
import json
from pathlib import Path
from protspect.sklearn import SpectralTransformer

transformer = SpectralTransformer(
    embedding=embedding,
    spectrum_length=8,
    standardize_channels=True,
)
X = transformer.fit_transform(["ACDE", "ACDF"])
provenance = {
    "aaindex": database.metadata(ids),
    "spectral_features": transformer.get_feature_metadata(),
}
Path(".aaindex/model-provenance.json").write_text(
    json.dumps(provenance, indent=2, allow_nan=False), encoding="utf-8"
)
```

The source metadata includes selected IDs/order, full original values, literature
citations, comments/correlation text, path, URL when supplied or downloaded,
version label, and raw file hash. Export this beside fitted preprocessing
metadata: transforming an embedding into ordinary dictionaries does not attach
AAindex's citation fields to the estimator automatically. Fit standardization
and any outcome-based index selection only within training folds.

Many indices are highly correlated. Select channels deliberately and begin
without cross spectra: 566 channels imply 160,461 within-sequence channel pairs
before frequency bins. This feature enables representations; it does not
establish improved predictive performance.

## Database attribution

Data comes from the [official AAindex1 download](https://www.genome.jp/ftp/db/community/aaindex/aaindex1).
See the [AAindex help and format specification](https://www.genome.jp/aaindex/aaindex_help.html).
Cite AAindex and the original publications of selected indices. One database
reference is Kawashima et al. (2008), *AAindex: amino acid index database,
progress report 2008*, Nucleic Acids Research 36, D202–D205,
[PMID 17998252](https://pubmed.ncbi.nlm.nih.gov/17998252/).

Protspect distributes the parser, not the AAindex database. Downloaded data
retains its own attribution and applicable terms; the package's BSD license
does not assert a license over downloaded databases. Tests use original synthetic
records and do not need the website. The complete official file was also parsed
and checked locally during development.
