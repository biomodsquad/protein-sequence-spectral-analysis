# Release record

## Protspect 0.1.0a2 — 2026-10-27

- [Production PyPI](https://pypi.org/project/Protspect/0.1.0a2/)
- [GitHub prerelease](https://github.com/biomodsquad/protein-sequence-spectral-analysis/releases/tag/v0.1.0a2)
- [Successful TestPyPI rehearsal](https://github.com/biomodsquad/protein-sequence-spectral-analysis/actions/runs/40589649064)
- [Successful production upload and verification](https://github.com/biomodsquad/protein-sequence-spectral-analysis/actions/runs/40590396808)
- [Cross-platform pre-merge checks](https://github.com/biomodsquad/protein-sequence-spectral-analysis/actions/runs/40587609566)
- [Machine-readable artifact and environment evidence](../validation/release_0.1.0a2.json)

Includes Kidera factors, AAindex1/2/3 adapters, Fourier/Z/Mellin/wavelet
transforms, unified sequence and paired collections, fitted sklearn pipelines,
and three bundled tutorial notebooks. AAindex data remains an explicit user
download. Wavelets require the optional `wavelets` extra.

```bash
python -m pip install 'Protspect[ml,wavelets]==0.1.0a2'
```

The production wheel and source archive are byte-for-byte the rehearsed files.
Both were downloaded from production PyPI, hash-verified, and installed in fresh
core-only and optional-dependency environments outside the checkout. Optional
stages passed 175 tests with one expected skip per artifact. The complete local
suite passed 188 tests; strict documentation, 15 doctests, and notebook execution
also passed. Alpha APIs may change.

| Artifact | SHA256 |
| --- | --- |
| Wheel | `5fe6e7a7c3431016d912c1a24c85b316f56220202c9b75b11ff24302295452cb` |
| Source archive | `ab5ce90289560bdb91971dd4d9936e91e889df98844a2ba5c285aa74721d6c61` |

Candidate: `abb28a8bc30ea7bf350f3362836c81f103ffb84f`; GitHub release tag:
`v0.1.0a2`, targeting `e2501d845f057393bd074a7d531bfe6870ae244d`.


## Protspect 0.1.0a1 — 2026-10-27

- [Production PyPI](https://pypi.org/project/Protspect/0.1.0a1/)
- [GitHub prerelease and artifacts](https://github.com/biomodsquad/protein-sequence-spectral-analysis/releases/tag/v0.1.0a1)
- [Production workflow](https://github.com/biomodsquad/protein-sequence-spectral-analysis/actions/runs/40544496824): candidate, publish, verify, and release jobs all passed.
- [Release workflow PR #3](https://github.com/biomodsquad/protein-sequence-spectral-analysis/pull/3): all 22 package checks passed before merge.
- [TestPyPI rehearsal](TESTPYPI.md): same wheel and source archive, built at `5f18766`.

The production workflow required the successful rehearsal and identical runtime
package code/metadata on main before publishing. Both production downloads
matched the rehearsed SHA256 hashes, then passed fresh core-only and optional
wheel/source installation checks outside the source checkout, numerical tests,
examples, and plotting. Dependencies were resolved from production PyPI.
GitHub tag `v0.1.0a1` points to `7a60a99259fe8399d0c16016e417b9a3f2ad19a0`,
which adds the production workflow without changing runtime package code.

| Artifact | SHA256 |
| --- | --- |
| Wheel | `71d49afe8e14b51f1dc1aade429acfdddc7f74a6360d4e8109b0cd1db21bd7cc` |
| Source archive | `3aa3edbe4c0c31c62bb2081546a67968e22c7d0e0aade77f7bb38a2398d47fe2` |

Install the alpha with:

```sh
python -m pip install Protspect==0.1.0a1
```

Do not rerun publishing for this version: uploaded distributions are immutable.
For a future release, update version/citation/changelog, validate a new candidate,
rehearse it on TestPyPI, and update the production workflow's pinned rehearsal
identity, version, and tag. Numerical corrections require affected historical
features to be recomputed and downstream models retrained. Alpha APIs may change.
