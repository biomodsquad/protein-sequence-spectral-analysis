# Protspect pre-publication review

Reviewed 2026-10-27. Original upstream commit:
`b31d4da5fbea0fdcc3c948429ddade0a7dfba3d2`.
Work is on local branch `codex/prospect-release-readiness`.

## Identity and scope

The GitHub repository was renamed from `biomodsquad/fourier-sequence-analysis`
to `biomodsquad/protein-sequence-spectral-analysis`; the local origin uses the
new URL. This checkout is independent of the surrounding language-model project.

PyPI's `prospect` name is already registered to Jacob Deppen's archaeological
survey simulation package. PyPI names are case-insensitive, so `ProSpect`
cannot be a separate distribution. Following the user's proposed replacement,
the draft distribution is **Protspect** and its import is **protspect**.
The PyPI JSON endpoint for `protspect` returned HTTP 404 on the review date.
This does not reserve the name or guarantee upload eligibility. No package was
uploaded and no changes were pushed or merged.

## Findings and implemented corrections

| Priority | Original behavior | Correction |
| --- | --- | --- |
| Critical | `compute_spectra` attached a centered frequency grid to unshifted FFT output. At length 4, a constant signal's DC coefficient appeared under -0.5 rather than 0; selected features could contain the wrong coefficients. | Center coefficients with `fftshift`, derive bins with `fftfreq`, and retain an explicitly documented periodic closure sample. Test direct DFT agreement for odd/even lengths, DC, Nyquist, phase, and Parseval's identity. |
| High | BLOSUM100 P-to-A was +2 while A-to-P was -2. | Correct P-to-A to -2; all 400 canonical entries were compared against the NCBI reference matrix and match. Test symmetry, expected dimensions, corrected entry, and immutability. |
| High | `zip` silently truncated mismatched chain collections and metadata columns. | Reject inconsistent row counts and reserved metadata names. Materialize generators once. |
| High | No build metadata, dependency declarations, release checks, or tests existed. | Add `pyproject.toml`, BSD license metadata, an alpha version, wheel/sdist packaging, numerical regression tests, and CI. |
| Medium | Empty sequences, unknown symbols, invalid lengths, and incorrect method order failed inconsistently. | Add explicit validation with descriptive errors. Preserve case-sensitive input rather than silently changing sequences. |
| Medium | Re-embedding or changing FFT length left stale derived spectra available. | Invalidate dependent state when inputs change; build sequence and paired collections before assigning their state. |
| Medium | `compute_cross_cross_spectra=True` was accepted but ignored. | Raise NotImplementedError for the unsupported operation. |
| Medium | Public collection classes and utility functions lacked documentation; the README contained no usable workflow. | Add module/class/method docstrings, runnable examples, numerical conventions, input policy, installation guidance, and migration notes. |
| Low | Misspelled positive-frequency selector and legacy package name made API discovery difficult. | Add correctly spelled selector and `protspect` namespace, retaining `FFTpro` imports and the misspelled alias. |

**Compatibility:** retaining old import paths does not preserve old feature
values. Recompute features and retrain models. There is no automatic conversion
of previously stored feature tables. The initial correction audit preserved raw
scaling. Later opt-in preprocessing options are documented below; raw defaults
remain unchanged.

## Validation

The subsequent ABP_GN audit adds four regression tests (28 total, passing),
16 fitted final classifiers across eight feature configurations, and paired
bootstrap comparisons. No convergence warnings were recorded. Its results are
mixed rather than a universal predictive improvement; see the scientific report.

- 24 tests pass on local Python 3.12 using NumPy 1.26.4, SciPy 1.13.1,
  and pandas 2.3.3.
- Ruff's configured correctness checks pass for package, tests, and the new
  self-contained example. Historical notebooks are excluded from this check.
- Both wheel and source distribution build with an isolated setuptools backend;
  Twine metadata/readme validation passes. The wheel contains the two package
  namespaces and license, excluding research notebooks and datasets.
- The renamed wheel was installed into a local `.venv` using the existing
  scientific dependencies. All 24 tests and the basic example also pass outside
  the checkout against that installed wheel; this was not a fresh dependency
  resolution or a TestPyPI install.
- The current 21-job installed-distribution CI matrix and fresh local minimum/
  latest dependency checks are documented in [COMPATIBILITY.md](docs/COMPATIBILITY.md).

## Normalization and grid implementation validation

The recommendation-2 implementation passes 75 tests, including direct DFT and
Parseval scaling identities, correctly scaled cross/paired spectra, zero-energy
and constant channels, fixed schema across batch lengths, training-only channel
statistics during cross-validation, sklearn cloning/GridSearchCV/Pipeline,
pandas output, and pickle/metadata round trips. The original defaults continue
to reproduce the pre-normalization outputs. The full ABP_GN dataset was also
transformed with training-fitted preprocessing: (3,562, 220) training features and
(1,973, 220) test features, all finite. Fitted metadata stayed unchanged after
test transformation. This is an integration check, not normalization selection.

## Recommended before publication

1. **Validate scientific results after the corrections.** An initial controlled
   ABP_GN audit is now available in [validation/ABP_GN_RESULTS.md](validation/ABP_GN_RESULTS.md),
   with all four correction combinations, original/default amplitude feature
   types, and the previous classifier-selection protocol. Follow-up should
   include repeated seeds and a new or family-disjoint test set. Rerun representative
   downstream analyses and compare predictive performance against the original
   features. The numerical tests establish transform correctness, not biological
   usefulness. Record embedding, FFT length, selector, version, and feature
   column order with exported datasets and trained models.
2. **Normalization and grid policies: implemented, scientific selection pending.**
   Opt-in length, square-root-length, and unit-energy scaling and mean centering
   preserve raw defaults. The optional sklearn SpectralTransformer learns
   channel statistics from training residues only and freezes its grid,
   selector, embedding, and output schema. It exports exact fitted metadata.
   See [docs/PREPROCESSING.md](docs/PREPROCESSING.md). The full AMP integration
   check preserves train/test columns without fitting on test data. Select any
   preferred normalization using training-only validation before changing a
   scientific workflow; this implementation claims no universally best policy.
3. **Research example: implemented using the complete official website data.**
   The archived notebook's outputs are cleared. A replacement source notebook and
   script use a pinned MISTIC 0.2.0b1 commit and an isolated, resolved Python 3.12
   environment. Website preparation verifies the pinned official download hashes,
   preserves every record and official identifier, and records raw/prepared-file
   provenance. The active example ignores the historical CSVs. Training contains
   1,146 records; official test contains 2,341, with 2,337 evaluated after excluding
   four exact positive training overlaps. See [research/README.md](research/README.md)
   and [results](research/RESULTS.md). Both source studies are cited. Raw downloads
   and prepared CSVs are local cache files excluded from distributions; research
   dependencies remain outside runtime requirements.
4. **Confirm the supported environment range — implemented and verified.**
   All 21 installed-distribution jobs passed, covering Python 3.10–3.15 on Linux,
   macOS, and Windows, with exact minimum stacks on Python 3.10. Fresh local
   installs found and corrected a scikit-learn lower-bound incompatibility.
   See [compatibility evidence](docs/COMPATIBILITY.md).
5. **Finalize release metadata and ownership — implemented.** The maintainer
   approved group authorship/maintenance, GitHub Issues as the public contact,
   a software citation, and `0.1.0a1` for the first alpha. These are recorded in
   package metadata, `CITATION.cff`, and README. Added changelog, support policy,
   and contribution guidance, all included in the source distribution. PyPI's
   JSON endpoint for `protspect` returned HTTP 404 on 2026-10-27; this does not
   reserve the name or guarantee registration. Recheck immediately before
   upload. A dedicated paired-chain result type remains a future API change.
6. **Perform a TestPyPI rehearsal — passed.** Trusted publishing uploaded
   `0.1.0a1` to TestPyPI. Downloaded wheel and source SHA256 hashes match the
   validated candidate; clean installations, tests, examples, and plotting
   passed with dependencies from production PyPI. Published metadata and local
   rendering of the published README passed inspection. TestPyPI's client
   challenge prevented visual page inspection; check the linked page manually.
   See [rehearsal evidence](docs/TESTPYPI.md). Production publishing remains a
   separate release decision.

## References

- [NCBI BLOSUM100 reference matrix](https://ftp.ncbi.nlm.nih.gov/blast/matrices/BLOSUM100).
- [Existing prospect project on PyPI](https://pypi.org/project/prospect/).
- [PyPI name normalization](https://packaging.python.org/en/latest/specifications/name-normalization/).
- [Writing pyproject.toml](https://packaging.python.org/en/latest/guides/writing-pyproject-toml/).
- [Packaging Python projects](https://packaging.python.org/en/latest/tutorials/packaging-projects/).
- [PyPI trusted publishing](https://docs.pypi.org/trusted-publishers/).


## Published alpha

`0.1.0a1` was published to production PyPI on 2026-10-27 after all release
checks passed. See [the production release record](docs/RELEASES.md). Earlier
review sections describe the work and decisions leading to that release.
