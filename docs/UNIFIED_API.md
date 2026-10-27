# Unified sequence-transform API (unreleased)

Use `method='fourier'|'z'|'mellin'|'wavelet'` throughout the primary API.
Every entry point uses the same embedding, numerical backend, coefficient-product
rules, and feature exporter. The PyPI 0.1.0a1 release predates this API; install
this checkout. Wavelets require `.[wavelets]`; sklearn requires `.[ml]`.

| Work | Primary interface |
| --- | --- |
| One sequence | `transform_sequence` or `Sequence` |
| Many sequences | `SequenceSet` |
| One two-chain sample | `transform_pair` |
| Many two-chain samples | `ComplexSet` |
| Fitted sklearn pipeline, single or paired samples | `SequenceTransformer` |

## Individual sequences

```python
from protspect import Sequence, transform_sequence
from protspect.embeddings import b100

result = transform_sequence('ACDE', b100, method='fourier', length=8)
values, names = result.features(('amplitude', 'phase'))

sequence = Sequence('ACDE')
sequence.embed_sequence(b100)
sequence.compute_transform(8, method='fourier')
values, names = sequence.features(('amplitude', 'phase'))
```

`TransformResult` contains `coefficients`, `coordinates`, and `metadata`.
`Sequence.result` is that same result type. Fourier/Z coordinates include
`frequency` in cycles per residue; Z also includes complex `z`. Mellin uses
`log_frequency`. Wavelets use `band` and within-band `index`. See
[definitions, defaults, and limitations](TRANSFORMS.md).

## SequenceSet remains a primary API

```python
from protspect import SequenceSet

sequences = SequenceSet(method='wavelet', compute_cross_spectra=False,
    transform_options={'wavelet': 'haar', 'level': 2, 'mode': 'periodization'})
sequences.from_sequence_list(['ACDE', 'ACDF'], length=8,
                             attributes={'id': ['p1', 'p2']})
table = sequences.export_features('real', attributes=['id'])
```

Use the same `method` for all four backends. `embedding`, `normalization`, and
`center` configure the collection. `transform_options` contains backend-specific
arguments (`radius`, `mellin_c`, `wavelet`, `level`, `mode`). The shared `length`
is explicit or the next power of two covering all chains. `results` provides the
computed TransformResults; `monomers` retains sequence objects and metadata.

`compute_cross_spectra` controls historical eager Fourier-array storage only.
Set it False to avoid this extra storage when using primary exports. New
`cross_coefficients` exports compute products on demand for every backend.

## ComplexSet remains a primary API

```python
from protspect import ComplexSet, transform_pair

pairs = ComplexSet(method='z', transform_options={'radius': 1.05})
pairs.from_sequence_list(['ACDE'], ['ACDF'], length=8,
                          attributes={'pair_id': ['complex-1']})
paired_table = pairs.export_features(('amplitude', 'phase'), attributes=['pair_id'])
chain_a_table = pairs.export_features('amplitude', chain_ID='A')
chain_b_table = pairs.export_features('amplitude', chain_ID='B')
result = transform_pair('ACDE', 'ACDF', method='z', length=8, radius=1.05)
```

A/B collections must have equal counts. Chains share one grid and are preprocessed
independently. `A`, `B`, and `AB` hold chain/pair objects; `results` holds paired
TransformResults. Every pair computes all ordered `A_i * conjugate(B_j)` channel
products. Reversing chains conjugates the corresponding reversed-channel product.
These are aligned coefficient descriptors. Wavelet products describe matched
band/position coefficients, not frequency-domain cross-spectral density or a
physical contact model. `compute_cross_cross_spectra=True` remains unsupported.

## Shared feature controls

`features` on a result/sequence and `export_features` on a set accept:

- `component`: amplitude, phase, real, imaginary, or an ordered tuple/list.
- `kinds`: `('coefficients',)`, `('cross_coefficients',)`, or both. Within-result
  products use the upper channel triangle, including auto-products. For an AB
  result its `coefficients` are already paired products; requesting another
  cross-product layer forms higher-order products.
- `selector`: a boolean mask or a callable receiving a **coordinate dictionary**.
  The mask must match coefficient length. Default None selects all coefficients.

```python
from protspect.utility import frequency_select_physical

# Actual Fourier/Z axes only; do not apply this to log-frequency or wavelet axes.
values, names = result.features(selector=lambda c: frequency_select_physical(c['frequency']))
# For a wavelet result, use selector=lambda c: c['band'] == 'D1'.
```

Names are channel-major, then component, then original coefficient index.
Cross-coefficient names include `cross_` before the component. Array indices
remain original indices after selection; inspect coordinates for their meaning.
Real wavelet coefficients retain signs; magnitudes lose them. Phase at zero
magnitude has no useful interpretation. Products inherit channel scaling and can
overflow if the embedding is too large; nonfinite results are rejected.

## One fitted sklearn implementation

```python
from protspect.sklearn import SequenceTransformer
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression

features = SequenceTransformer(method='fourier', length=32,
    component=('amplitude', 'phase'), kinds=('coefficients', 'cross_coefficients'),
    selector=lambda c: c['frequency'] >= 0,
    normalization='length', standardize_channels=True)
model = make_pipeline(features, StandardScaler(), LogisticRegression())
# model.fit(training_sequences, labels)

paired_features = SequenceTransformer(method='wavelet', length=32, paired=True,
    component='real', transform_options={'wavelet': 'haar', 'level': 3})
X = paired_features.fit_transform([('ACDE', 'ACDF'), ('AAAA', 'CCCC')])
metadata = paired_features.get_feature_metadata()
```

Fit inside each CV fold. Length and channel moments are learned from training
chains only; with paired=True both training chains contribute, one weight per
residue. Selection is evaluated once at fit, frozen, and reused on later batches.
The fitted grid, components, options, embedding, and column order are snapshots.
Longer prediction chains raise instead of expanding the schema. Changing parameters
with `set_params` invalidates fit; direct constructor mutations require refitting.
Metadata is independent and JSON-safe, including real/imaginary parts of Z axes.
Labels are ignored by the featurizer. `.set_output(transform='pandas')` is supported.

## Compatibility and migration

`SequenceFFT` is a compatibility subclass of `Sequence`; `SpectralTransformer`
is a constructor/metadata adapter over `SequenceTransformer`. They share the
numerical/export/fitting implementation, rather than parallel FFT code paths.
`FFTpro` imports continue working. Fitted snapshots from the preceding adapters
are upgraded on unpickle; no deprecation warnings are emitted yet.

| Existing call | Primary equivalent |
| --- | --- |
| `SequenceFFT(...).compute_spectra(L)` | `Sequence(...).compute_transform(L, method='fourier')` after embedding |
| `compute_spectral_properties(...)` | `features(component, kinds=..., selector=...)` |
| `from_sequence_list(..., spectrum_length=L)` | `from_sequence_list(..., length=L)` |
| `export_datatable(...)` | `export_features(...)` |
| `SpectralTransformer(spectrum_length=L, ...)` | `SequenceTransformer(method='fourier', length=L, ...)` |
| `transform_sequence(..., transform='wavelet')` | `transform_sequence(..., method='wavelet')` |

Existing Fourier calls preserve their original values, Amp_/Phase_ names,
frequency selectors, metadata schema, and L+1 periodic closure. Primary exports
use L unique Fourier bins, indexed names, and **all bins by default**; they do
not silently inherit the old dyadic selector. For even lengths the primary grid
contains -0.5, not the repeated +0.5 endpoint. Thus old and primary feature tables
are not interchangeable by column position: retain compatibility calls for an
existing model, or recompute and refit with the primary API. This migration does
not change the AMP notebook's full-coefficient feature definition.

Old `export_datatable` remains a Fourier-only compatibility method and raises
for other backends with guidance to use `export_features`. The set classes
themselves are shared across all backends. `spectrum_length` stays as a spelling
alias for collection input length; supplying both spellings raises. Earlier AMP
baseline reports describe their original source snapshots; reruns have their own
source hashes and output folders.
