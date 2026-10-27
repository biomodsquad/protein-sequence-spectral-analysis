# Protspect

![Protspect logo](docs/_static/protspect-logo.png)

**Protein sequence spectral analysis**: transform features for proteins and
peptides using numeric residue-property channels. This is the original sequence
featurization method, separate from the subsequent protein language-model work.

The repository is `biomodsquad/protein-sequence-spectral-analysis`. The
[PyPI distribution](https://pypi.org/project/Protspect/) is **Protspect**, imported as **`protspect`**. The original
`ProSpect` distribution name conflicts with an unrelated existing PyPI package.
Version `0.1.0a1` is the first published alpha; APIs may change.

## Documentation and tutorials

The [documentation source](docs/index.rst) follows the MISTIC Sphinx/Furo layout:
installation, matrix-based concepts, eight guided tutorials, downloadable
notebooks, and a public API reference. The supplied logo and matrix illustration
introduce the residue embedding → zero-padding → sequence-axis transform workflow.

Three self-contained notebooks are bundled in `protspect.examples`; locate them
with `importlib.resources.files("protspect.examples")` and copy to a writable
workspace. They cover matrix-to-feature calculations, four transforms and paired
chains, and training-only sklearn pipelines. The full AMP MLP comparison and
embedding atlas remain separate research examples with their documented inputs.

Build a local site with:

```sh
python -m pip install -e '.[docs]'
python -m docs.generate_tutorial_figures
python -m sphinx -W --keep-going -b html docs docs/_build/html
```

See [build and validation instructions](docs/contributing.rst) for executable
snippets, notebook checks, and manual Pages publication. The documentation
explicitly separates the development API from the published 0.1.0a1 alpha.

## Install

Requires Python 3.10 or later. Install the published alpha and its dependencies:

```bash
python -m pip install Protspect==0.1.0a1
# Optional scikit-learn integration:
python -m pip install 'Protspect[ml]==0.1.0a1'
```

For development, install from the repository:

```bash
git clone https://github.com/biomodsquad/protein-sequence-spectral-analysis.git
cd protein-sequence-spectral-analysis
python -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[dev]'
```

The distribution is named `Protspect`. The unrelated `prospect` package is a
different project.

## Alternative transforms (development checkout)

Fourier, circular Z, finite discrete Mellin, and discrete wavelet transforms now
share `transform_sequence`, with a fixed-width `SequenceTransformer` for sklearn.
See [definitions and examples](docs/TRANSFORMS.md). Wavelets require the optional
`wavelets` extra. These additions are not yet in the published alpha.

## Residue representations

The development branch includes the ten Kidera factors (KF1–KF10), in addition
to the default BLOSUM100 embedding. See [representation documentation and
citation](docs/REPRESENTATIONS.md). This addition is not in PyPI version 0.1.0a1.

```python
from protspect import SequenceSet
from protspect.embeddings import kidera_factors

features = SequenceSet(embedding=kidera_factors(), compute_cross_spectra=False)
features.from_sequence_list(["ACDE", "ACDF"], spectrum_length=8)
```

AAindex1 files can also supply selected property channels (development version):

```python
from protspect.aaindex import AAIndex

database = AAIndex.from_file(".aaindex/aaindex1")
features = SequenceSet(
    embedding=database.embedding(["KYTJ820101", "FAUJ880103"]),
    compute_cross_spectra=False,
)
```

See [AAindex1 loading, explicit downloads, and provenance](docs/AAINDEX.md).
No database is bundled or downloaded during feature extraction.

AAindex2/3 matrices also support explicit row or column profile channels:

```python
from protspect.aaindex import AAIndex2

matrices = AAIndex2.from_file(".aaindex/aaindex2")
embedding = matrices.embedding("HENS920102", axis="rows")  # BLOSUM62
```

See [matrix conversion and contact-potential interpretation](docs/AAINDEX_MATRICES.md).

## One transform API

The development API supports `method='fourier'|'z'|'mellin'|'wavelet'` throughout.
`SequenceSet` and `ComplexSet` remain primary interfaces for single and paired
collections. All entry points share the same backend and feature exporter.
These additions are not yet in PyPI 0.1.0a1; install this checkout, adding
`.[wavelets]` for wavelets and `.[ml]` for sklearn.

```python
from protspect import SequenceSet, ComplexSet, transform_sequence

features = SequenceSet(method='fourier', compute_cross_spectra=False)
features.from_sequence_list(['ACDE', 'ACDF'], length=8,
    attributes={'name': ['peptide-1', 'peptide-2']})
table = features.export_features('amplitude', attributes=['name'])
print(table.shape)  # (2, 161): metadata + 20 channels × 8 unique Fourier bins

pairs = ComplexSet(method='z', transform_options={'radius': 1.05})
pairs.from_sequence_list(['ACDE'], ['ACDF'], length=8)
paired_table = pairs.export_features(('amplitude', 'phase'))
chain_a_table = pairs.export_features('amplitude', chain_ID='A')

result = transform_sequence('ACDE', method='fourier', length=8)
values, names = result.features(selector=lambda c: c['frequency'] >= 0)
```

The default BLOSUM100 embedding supplies 20 residue-property channels. Kidera,
AAindex, and custom channel mappings work with every transform. Inputs are
case-sensitive; unknown residues are rejected when not covered by the embedding.
All chains in a collection share a fixed grid; no sequence truncation occurs.
Per-sequence centering and none/length/sqrt_length/energy normalization operate
on real residues before padding. Metadata cannot overwrite features.

`component` selects amplitude, phase, real, imaginary, or an ordered tuple;
`kinds` selects coefficients and/or within-channel `cross_coefficients`.
`selector` receives named coordinate arrays: frequency for Fourier/Z,
log-frequency for Mellin, band/index for wavelets. Primary exports include
all coefficients by default, with no periodic duplicate. Paired products are
aligned coefficient descriptors, not a physical model of residue contacts.

For fitted pipelines use `protspect.sklearn.SequenceTransformer`, with the same
method/options/components and coordinate selectors. Set `paired=True` for
(A, B) samples. Grid and channel moments fit on training chains only; selection
and output columns are frozen. See [the unified API](docs/UNIFIED_API.md),
[transform definitions](docs/TRANSFORMS.md), and
[normalization formulas](docs/PREPROCESSING.md).

`SequenceFFT`, `SpectralTransformer`, and the sets' `export_datatable` methods
remain compatibility adapters. Existing Fourier values, selectors, closure bins,
and feature names are preserved. Use them with existing fitted models; recompute
and refit before switching to the primary full-bin feature schema.

## Migration from FFTpro

`FFTpro.features`, `FFTpro.embeddings`, and `FFTpro.utility` remain compatibility
imports. New code should use `protspect`. The misspelled selector
`frequency_selelect_positive` is retained as an alias.

**Feature values have intentionally changed.** The original implementation
labeled unshifted FFT coefficients using a centered grid, assigning coefficients
to incorrect frequencies. The corrected implementation aligns bins and labels.
The BLOSUM100 P-to-A score was also corrected from +2 to -2, matching A-to-P.
Regenerate existing feature tables and retrain downstream models; old and new
features must not be mixed. Empty sequences, short FFT lengths, unequal pair
counts, and unequal metadata lengths now raise errors instead of failing later
or silently dropping records.

The historical notebook in `FFTpro/examples/` is archived with all saved outputs
cleared; the original CPP CSVs remain unchanged and excluded from distributions.
A new [reproducible CPP notebook and workflow](research/README.md) uses the
complete official MLCPP 2.0 website datasets and a pinned MISTIC source, a separate locked research environment, training-only model
selection, exact train/test overlap checks, and recorded provenance. See
[executed results](research/RESULTS.md) and
[official data preparation and provenance](research/DATA_PROVENANCE.md).
Research dependencies are optional and separate from the runtime package.
Use `examples/basic_features.py` for a minimal core-only example.

## Development and release

```bash
python -m pytest -q
python -m ruff check .
python -m build
python -m twine check dist/*
```

For the controlled antimicrobial-peptide correction audit, see
[validation/README.md](validation/README.md). Install `.[dev,validation]` to run
its classifier tests and benchmark.

CI checks Python 3.10–3.15 on Linux, macOS, and Windows, installing both
wheel and source archives into fresh environments outside the checkout. Exact
minimum dependencies are checked on Python 3.10. See
[the compatibility policy](docs/COMPATIBILITY.md) for scope and evidence. See [RELEASE_REVIEW.md](RELEASE_REVIEW.md) for findings,
validation evidence, and the remaining release decisions. The initial alpha version is
`0.1.0a1`. The [release record](docs/RELEASES.md) links the successful TestPyPI rehearsal
and production download/install checks.

## Citation and project contact

Cite the software version and commit used:

> Biomolecular Systems Modeling and Design Group. *Protspect: protein sequence
> spectral analysis*. Version 0.1.0a1.
> https://github.com/biomodsquad/protein-sequence-spectral-analysis

Machine-readable metadata is in [CITATION.cff](CITATION.cff). Cite dataset
sources and external methods separately; the CPP example lists these in
[its provenance documentation](research/DATA_PROVENANCE.md).

The group maintains this package through [GitHub Issues](https://github.com/biomodsquad/protein-sequence-spectral-analysis/issues).
See the [changelog](CHANGELOG.md), [support policy](SUPPORT.md), and
[contribution guide](CONTRIBUTING.md).

Licensed under the [BSD 3-Clause License](LICENSE).
