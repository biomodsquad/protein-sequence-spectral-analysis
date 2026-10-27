# Changelog

## 0.1.0a2 — 2026-10-27

- Add a Sphinx/Furo documentation site with eight guided tutorials, generated
  figures, API reference, three bundled executable notebooks, and docs CI.
  Include the supplied project logo and matrix/sequence-transform illustration.

- Unify single sequences, SequenceSet, ComplexSet, and sklearn around one
  transform/export engine. Add coordinate selectors, cross-coefficient exports,
  paired pipelines, and shared metadata. Preserve old Fourier APIs as adapters.

- Add explicit Fourier/Z/Mellin/discrete-wavelet coefficient transforms, optional
  PyWavelets dependency, and fixed-schema SequenceTransformer integration.

- Add AAindex2/3 matrix catalogues, explicit section download helpers, native
  axis/storage provenance, and canonical row/column profile adapters. Preserve
  asymmetric/rectangular matrices and reject missing canonical scores.


- Add offline AAindex1 parsing, search, channel selection, explicit official
  downloads, SHA256 pinning, and selected-record provenance. Missing values
  are rejected during embedding conversion; no database is redistributed.

- Add all ten Kidera residue-factor channels, with `kidera_factors()` for
  independent embedding dictionaries and `kidera` for read-only access.
  BLOSUM100 remains the default. See [representation provenance](docs/REPRESENTATIONS.md).

## 0.1.0a1 — 2026-10-27

First packaged alpha of Protspect, the protein/peptide spectral featurization
method previously developed as FFTpro. Published to PyPI after the TestPyPI
rehearsal; see [the release record](docs/RELEASES.md).

### Added

- Installable `protspect` package with optional ML and validation dependencies.
- Public API docstrings, numerical regression tests, and feature examples.
- Length, square-root-length, and energy normalization, sequence centering,
  and a fitted scikit-learn transformer with fixed feature grids, training-only
  channel statistics, named outputs, and exported preprocessing metadata.
- Reproducible CPP analysis using official MLCPP2 downloads and pinned MISTIC
  source, with split provenance and separate validation/test evaluation.
- Fresh wheel/source installation checks on CPython 3.10–3.15 across Linux,
  macOS, and Windows, including exact minimum dependency stacks on 3.10.

### Corrected

- Align FFT frequency bins with shifted coefficients.
- Correct the BLOSUM100 proline-to-alanine entry from +2 to -2.
- Validate inputs and paired-chain counts; invalidate stale features when
  sequences change; reject the unimplemented cross-cross spectrum option.

### Migration and compatibility

`FFTpro` imports remain compatibility wrappers. Raw defaults are preserved,
but corrected frequency bins and embedding scores change affected feature
values. Recompute affected feature tables and retrain downstream models; do
not mix historical and corrected features. Record the package version and
preprocessing metadata with trained models.

The optional scikit-learn minimum is 1.4 because older output wrapping fails
for list input with pandas output. Alpha APIs may change between releases.
A dedicated paired-chain result type is deferred to a later API revision.
