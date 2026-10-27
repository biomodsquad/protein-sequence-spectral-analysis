# TestPyPI rehearsal

This workflow publishes only to TestPyPI. It does not publish to production
PyPI, create a GitHub release, merge the review PR, or tag a release.

## Account-side setup

Sign in to [TestPyPI publishing](https://test.pypi.org/manage/account/publishing/)
and create a pending GitHub publisher for this new project:

| Field | Value |
| --- | --- |
| Project name | `Protspect` |
| GitHub owner | `biomodsquad` |
| Repository | `protein-sequence-spectral-analysis` |
| Workflow filename | `release-testpypi.yml` |
| Environment | `testpypi` |

TestPyPI uses a separate account and publisher configuration from production
PyPI. Do not store an API token in the repository. The GitHub `testpypi`
environment permits the review branch and `main`; other branches cannot deploy.
Only the publishing job has `id-token: write`, and it runs no project code.

## Workflow and checks

The workflow runs when first added/updated on the review branch and can be
rerun from Actions. Once present on the default branch it can also be dispatched
manually. Do not rerun the publishing job after a successful upload: release
files cannot be replaced. Rerun only verification if an upload already succeeded.

1. Build the candidate, check metadata, run package tests, and test fresh wheel
   and source installations before allowing publication.
2. Upload exactly that build artifact to `https://test.pypi.org/legacy/` through
   OIDC trusted publishing. The destination is fixed in the workflow.
3. Download both artifacts from TestPyPI's version-specific JSON API, verify
   their SHA256 digests against the original build, and retain the server metadata.
4. Install each downloaded artifact in a separate fresh virtual environment
   outside the checkout. Dependencies come only from production PyPI; no mixed
   package indexes are used. Check core imports, numerical tests, optional
   dependencies, the full test suite, examples, and a plotting smoke test.
5. Inspect the rendered [TestPyPI project page](https://test.pypi.org/project/Protspect/0.1.0a1/)
   and its maintainer, license, Python/dependency bounds, project links, and README.

`testpypi-rehearsal-evidence` retains the downloaded distributions, JSON
metadata, exact environments, and test logs. Record the successful run URL and
commit in this document after completion.

## Status — rehearsal passed

[Rehearsal run 40537066711, attempt 3](https://github.com/biomodsquad/protein-sequence-spectral-analysis/actions/runs/40537066711/attempts/3)
passed candidate, publish, and verify jobs on 2026-10-27. The published
candidate was built from commit `5f18766`; later commits document setup/results.

- TestPyPI: [Protspect 0.1.0a1](https://test.pypi.org/project/Protspect/0.1.0a1/).
- Candidate job: 75 checkout tests, metadata validation, and clean wheel/source
  installation checks passed.
- Post-upload: SHA256-verified TestPyPI downloads matched the candidate exactly.
  Fresh wheel and source installs each passed 74 optional-suite tests with one
  expected Git-history skip, plus core-only tests, both examples, and plotting.
  Dependencies were resolved only from production PyPI on CPython 3.12.15.
- Published JSON metadata confirms name/version, group author/maintainer,
  BSD-3-Clause license, Python >=3.10, dependency bounds, and project URLs.
  The published README renders locally using `readme_renderer[md]`. TestPyPI's
  page returned a client challenge to inspection tools, so visual page layout
  has not been verified. Inspect the linked page manually.
- No production PyPI upload, merge, tag, or GitHub release was performed.

| Artifact | SHA256 |
| --- | --- |
| Wheel | `71d49afe8e14b51f1dc1aade429acfdddc7f74a6360d4e8109b0cd1db21bd7cc` |
| Source archive | `3aa3edbe4c0c31c62bb2081546a67968e22c7d0e0aade77f7bb38a2398d47fe2` |

Earlier attempts stopped before upload: first no matching publisher was
registered; then the registered project name differed from `Protspect`. After
correcting the account-side registration, attempt 3 succeeded. Do not rerun
the publishing job for this version; uploaded files are immutable.


This section records the initial `0.1.0a1` rehearsal. The subsequent `0.1.0a2`
rehearsal and production verification both passed; see [release records](RELEASES.md).
Do not rerun uploads for either existing version: uploaded files are immutable.
