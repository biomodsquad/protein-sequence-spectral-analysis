# Compatibility policy and validation

Protspect requires Python >=3.10. The CI target is CPython 3.10–3.15 on
Linux, macOS, and Windows. A support target is not a claim of completed
verification; the pull request checks record the actual platform results.

Each of the 18 latest-compatible dependency jobs installs both the wheel and
source distribution in separate fresh virtual environments outside the checkout.
It checks core-only imports, the legacy namespace, numerical tests, examples,
then optional ML/validation dependencies, the full test suite, and plotting.
Three additional Python 3.10 jobs use exact minimum versions on each OS.
Logs and resolved dependency versions are uploaded as CI artifacts. One test
requiring original Git history skips outside a checkout; it runs in the build job.

## Dependency bounds

| Dependency | Minimum | Reason/evidence |
| --- | --- | --- |
| NumPy | 1.23 | Exact minimum tested with Python 3.10 |
| SciPy | 1.9 | Exact minimum tested with Python 3.10 |
| pandas | 1.5 | Exact minimum tested with Python 3.10 |
| scikit-learn (optional) | 1.4 | 1.2 fails pandas output for list input; 1.4 passes |
| Matplotlib (validation) | 3.6 | Exact minimum plotting smoke test |
| threadpoolctl (validation) | 3 | Exact minimum tested |
| setuptools (build) | 77 | PEP 639 license metadata; first available release 77.0.1 tested |

The scikit-learn 1.2 failure is triggered by
`SpectralTransformer().set_output(transform="pandas")` with a list of sequences:
its output wrapper treats the list's `index` method as a pandas index. The
minimum is raised to 1.4, whose pandas adapter handles list input correctly.

No speculative upper bounds are imposed. CI resolves a coherent latest stack
for each interpreter; arbitrary mixtures of old and new dependency versions,
PyPy, free-threaded builds, and additional architectures are not established by
this matrix. The research notebook has its own pinned Python 3.12 environment.

## Local evidence

Fresh wheel and source installs passed on macOS arm64 with CPython 3.14.10
and 3.15.2 using current dependencies. The minimum Python 3.10.22 wheel check
passed 74 tests with one expected checkout-only skip. The source minimum check
is being completed after correcting an invalid setuptools 77.0.0 test pin.
Remote platform verification is pending the first pull request run.

Run the same checks after building with:

```sh
python ci/check_installed.py --output-dir compatibility-results
python ci/check_installed.py --constraints ci/constraints/minimum.txt --output-dir minimum-results
```
