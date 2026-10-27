# Contributing

Report reproducible problems or propose changes through
[GitHub Issues](https://github.com/biomodsquad/protein-sequence-spectral-analysis/issues).
For substantial changes, describe the scientific motivation and intended API
before implementation. Submit pull requests against `main`.

## Development

Use Python 3.10 or later in a virtual environment:

```sh
python -m pip install -e '.[dev,validation]'
python -m ruff check .
python -m pytest -q
python -m build
python -m twine check dist/*
```

The pinned research environment and notebook reproduction instructions are in
[research/README.md](research/README.md); research dependencies are optional
and are not runtime requirements.

For numerical changes, provide an independently derived expected result
(e.g. a direct DFT, known matrix value, or mathematical invariant). Explain
which existing features and trained models are affected. Fit preprocessing
only on training data and retain a held-out evaluation when reporting model
performance. Avoid tests that merely repeat the implementation.

Document public behavior and update [CHANGELOG.md](CHANGELOG.md). Preserve
legacy imports where practical; call out intentional compatibility changes.
Do not commit credentials, private sequences, downloaded dataset caches,
virtual environments, or recomputable feature matrices. Cite data sources and
verify redistribution terms before adding datasets.

Contributions to the package are distributed under its existing BSD-3-Clause
license. Keep third-party attribution and licensing separate from package
authorship; dataset and research dependencies retain their own terms.
