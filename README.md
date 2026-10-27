# Protspect

**Protein sequence spectral analysis**: Fourier features for proteins and
peptides using numeric residue-property channels. This is the original sequence
featurization method, separate from the subsequent protein language-model work.

The repository is `biomodsquad/protein-sequence-spectral-analysis`. The
[PyPI distribution](https://pypi.org/project/Protspect/) is **Protspect**, imported as **`protspect`**. The original
`ProSpect` distribution name conflicts with an unrelated existing PyPI package.
Version `0.1.0a1` is the first published alpha; APIs may change.

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

## Generate sequence features

```python
from protspect import SequenceSet

features = SequenceSet(compute_cross_spectra=False)
features.from_sequence_list(
    ["ACDE", "ACDF"],
    spectrum_length=8,
    attributes={"name": ["peptide-1", "peptide-2"]},
)
table = features.export_datatable(attributes=["name"])
print(table.shape)  # (2, 81): metadata plus 20 channels at four frequencies
```

The default embedding `protspect.embeddings.b100` maps 20 canonical amino acids
to their 20 BLOSUM100 substitution-score channels. Inputs are case-sensitive;
unknown residues, gaps, stop symbols, and lowercase letters are rejected when
not covered by the embedding. Custom numeric channels can encode other residue
properties or alphabets:

```python
from protspect import SequenceFFT

sequence = SequenceFFT("ACAC")
sequence.embed_sequence({"property": {"A": 1.0, "C": -1.0}})
sequence.compute_spectra(length=4)
properties = sequence.compute_spectral_properties(compute_phase=True)
assert properties["Amp_property_0.5"] == 4.0
```

`embed_sequence` stores one numeric vector per channel; `compute_spectra`
right-pads each vector with zeros and computes its FFT. `compute_cross_spectra`
optionally creates all upper-triangle within-sequence channel cross-products,
including auto-spectra. To export those, use
`features.export_datatable(spectra=("cross_spectra",))` after enabling them.
Metadata names must not collide with API attributes or generated feature names.

## Paired chains

```python
from protspect import ComplexSet

pairs = ComplexSet()
pairs.from_sequence_list(["ACDE"], ["ACDF"], spectrum_length=8)
cross_chain_features = pairs.export_datatable(chain_ID="AB")
chain_a_features = pairs.export_datatable(chain_ID="A")
```

A/B collections must contain the same number of sequences. Each pair produces
all ordered channel products `FFT(A_i) * conjugate(FFT(B_j))`. Swapping A and B
conjugates the corresponding product, so phase encodes chain order. These are
spectral descriptors, not a physical model of binding or residue contacts.
`compute_cross_cross_spectra=True` is explicitly unsupported.

## Numerical conventions

- All sequences in a collection share one FFT length. If omitted, the next
  power of two covering the longest sequence is used. Supply the same explicit
  length for training, validation, and test collections to preserve columns.
- FFTs are unnormalized by default and their bins are centered with `fftshift`. For FFT length `L`,
  `spectrum_length` is `L`, while each stored spectrum and its `frequency` array
  contain `L + 1` entries. The first bin is repeated at the end to close the
  periodic spectrum, with its frequency increased by one cycle per residue.
  For even `L`, this gives frequencies from -0.5 through +0.5. For odd `L`,
  the endpoint exceeds +0.5; it remains a periodic duplicate, not a new bin.
- The default selector exports DC and positive dyadic frequencies on the grid:
  1/2, 1/4, 1/8, etc. For non-power-of-two lengths some or all dyadic bins
  may be absent. A custom selector receives the frequency array and must return
  a boolean mask of the same shape. `frequency_select_positive` in
  `protspect.utility` selects all nonnegative entries, including the closure bin.
- `Amp_...` features are absolute magnitudes; `Phase_...` features are angles
  in radians. Cross-spectrum magnitudes scale quadratically with signal scale.
  Optional length, square-root-length, or unit-energy scaling and per-sequence
  mean centering are available before padding. No preprocessing is applied by
  default, and phase unwrapping is not implemented.
  Phase at zero magnitude is not meaningful.
- Padding, sequence length, and channel scale affect the representation.
  Amplitudes lose absolute circular-shift position information; phases depend
  on alignment and sequence origin. Standardize features using training data
  only when fitting downstream models.

## Normalization and fitted grids

`SequenceFFT.compute_spectra`, `SequenceSet`, and `ComplexSet` accept explicit
`normalization='none'|'length'|'sqrt_length'|'energy'` and `center=False|True`
options. Existing raw defaults are preserved.

For machine learning, install `.[ml]` and use
`from protspect.sklearn import SpectralTransformer` inside a scikit-learn
Pipeline. It freezes the grid, selector, embedding, and feature column order at
fit time, and can learn channel scales from training residues only. Longer
prediction sequences raise an error rather than changing the grid. See
[the preprocessing guide](docs/PREPROCESSING.md) for formulas, tradeoffs,
leakage prevention, and model/feature metadata persistence.

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
