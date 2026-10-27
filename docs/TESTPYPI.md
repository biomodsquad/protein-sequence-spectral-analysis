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

## Status

[Rehearsal run 40537066711](https://github.com/biomodsquad/protein-sequence-spectral-analysis/actions/runs/40537066711)
for commit `5f18766` passed the candidate job: build, Twine metadata checks,
75 checkout tests, and fresh core/optional wheel and source installation tests.
The local macOS CPython 3.12.2 rehearsal also passed both artifact checks.

TestPyPI rejected the OIDC exchange with `invalid-publisher`: no corresponding
publisher was found for the expected repository, workflow, and environment.
No upload occurred and post-upload verification was skipped. Configure the
account-side pending publisher above, then rerun failed jobs for this same run
to publish the already-tested candidate and perform post-upload verification.
Do not describe the end-to-end rehearsal as complete until those checks pass.

The current candidate is `0.1.0a1`, approved for the initial alpha. Before a
production release, recheck name availability and finish scientific/release
review. Production publishing is a separate action.
