# Release record

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
