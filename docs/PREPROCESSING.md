> Development API: use `SequenceTransformer` and the generic `SequenceSet` /
> `ComplexSet` exports for all four transforms. The Fourier examples below retain
> compatibility names and closed-bin conventions. See [the unified API](UNIFIED_API.md).

# Normalization and fitted feature grids

Raw, uncentered spectra remain the default. These options specify different
scientific representations; none is declared universally better. The ABP_GN
correction audit used raw spectra and does not select a normalization policy.

## Order of operations

1. Map actual residues to numeric channel scores.
2. Optionally standardize channels using **training-residue** means and population
   standard deviations, frozen by the fitted transformer.
3. Optionally subtract each channel's mean over the current sequence's residues.
4. Optionally scale each channel using one of the policies below.
5. Right-pad with zeros to the fixed FFT length, then apply the FFT and fftshift.
6. Form channel cross-products if requested, then export magnitudes and/or phases.
7. Optionally fit a downstream feature scaler within the training pipeline.

Padding is excluded from channel statistics, sequence means, and divisors.
Sequence centering and centering the FFT bin order with fftshift are distinct.
There is no windowing, detrending beyond mean subtraction, or phase unwrapping.

## Explicit normalization policies

For residue count `N` and channel vector `x` after optional centering:

| `normalization` | Signal entering padding | Interpretation |
| --- | --- | --- |
| `none` (default) | `x` | Preserve raw length and property-scale effects. |
| `length` | `x / N` | DC becomes the mean channel value when sequence centering is off. Does not eliminate every length effect. |
| `sqrt_length` | `x / sqrt(N)` | Signal energy becomes mean squared channel value; does not produce unit energy. |
| `energy` | `x / sqrt(sum(x**2))` | Unit-energy channel when nonzero; zero-energy channels remain exactly zero. |

`energy` divides by the **L2 norm**, not by squared energy and not by FFT
length. It uses a scaled calculation to avoid norm overflow/underflow.
Constant channels become exactly zero when sequence centering is enabled.

These policies are applied before forming cross-products. If channel i is
divided by `d_i` and j by `d_j`, `X_i * conjugate(X_j)` is divided by
`d_i * d_j`. Within one sequence, length normalization therefore divides
cross-spectra by `N**2`. For paired chains it divides by `N_A * N_B`.
The auto-spectrum is squared magnitude, not magnitude. No second normalization
is applied during export. Positive scaling preserves phases for nonzero bins;
centering may change them. Phase at zero magnitude is undefined scientifically,
even though numpy.angle supplies a numerical value.

```python
from protspect import SequenceFFT, SequenceSet, ComplexSet
from protspect.embeddings import b100

sequence = SequenceFFT('ACDE')
sequence.embed_sequence(b100)
sequence.compute_spectra(8, normalization='length', center=True)

# Same options are applied to every sequence or chain.
sequences = SequenceSet(normalization='energy', center=False)
sequences.from_sequence_list(['ACDE', 'ACDF'], spectrum_length=8)
pairs = ComplexSet(normalization='length')
pairs.from_sequence_list(['ACDE'], ['ACDF'], spectrum_length=8)
```

The low-level API has no fitted population statistics. It is appropriate when
using a predefined grid and fixed embedding. For machine-learning workflows,
use the fitted transformer below to preserve the entire feature schema.

## A fitted transformer for machine learning

Install `Protspect[ml]` (or `python -m pip install -e '.[ml]'` in a checkout), then:

```python
from protspect.sklearn import SpectralTransformer
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

model = make_pipeline(
    SpectralTransformer(
        spectrum_length=1024,
        normalization='length',
        center=False,
        standardize_channels=True,
        spectra_types=('spectra',),
    ),
    StandardScaler(),
    LogisticRegression(random_state=7),
)
model.fit(train_sequences, train_labels)
test_probabilities = model.predict_proba(test_sequences)[:, 1]
```

Only `protspect.sklearn` needs scikit-learn; importing the core `protspect`
package does not. The transformer consumes a one-dimensional iterable of strings
and returns a float64 matrix. To get a DataFrame directly, use
`SpectralTransformer(...).set_output(transform='pandas')`.
For a runnable small example, see `examples/normalized_pipeline.py`.

`fit` copies the embedding, chooses the FFT grid, learns optional channel
statistics, evaluates the frequency selector **once**, and freezes the ordered
output feature names. `transform` reuses these snapshots without learning from
the current batch. Consequently, the same sequence has the same representation
whether processed alone or alongside other sequences. Stateful selector changes
and later edits to the caller's embedding do not alter the fitted schema.

With `spectrum_length=None`, only the **training** sequences determine the next
power-of-two grid. A later sequence longer than this grid raises ValueError;
there is no silent truncation, interpolation, or grid expansion. To support
longer deployment inputs, set an explicit application-appropriate grid and refit
the entire pipeline. An input length rule selected before evaluation is preferable
to choosing it by inspecting held-out examples. `SequenceSet` still infers its
grid anew on each population call, so do not independently infer train/test grids.

Use `set_params(...)` and then `fit(...)` to change a fitted transformer. Calling
`set_params` invalidates its fitted state; a failed refit also invalidates the
old model. Direct constructor-parameter edits take effect on the next fit.
Paired-chain fitted transformation is not implemented; ComplexSet still requires
an explicit grid shared between training and prediction.

## Channel standardization and leakage

`standardize_channels=True` learns each channel's mean and population standard
deviation across **all residues in the fitting sequences**. Every residue has
equal weight, so longer sequences contribute more; padding contributes nothing.
A zero-variance training channel uses scale 1. The embedding still covers its
complete alphabet, including valid residues absent from the training sample;
unrecognized symbols raise errors. With standardization disabled, channel
means/scales are 0/1 and the original score mapping is preserved.

This is distinct from `StandardScaler` on the resulting feature columns. Keep
both inside the Pipeline when using GridSearchCV or cross-validation, so every
fold learns them only from its training partition. Passing precomputed globally
standardized features into cross-validation would leak validation information.
The transformer's labels argument is ignored.

Centering removes compositional means and DC, which may remove useful biological
signal. Unit-energy scaling removes channel energy and can amplify weak signals.
Combining channel standardization, sequence centering, and unit-energy scaling
can make some scaling operations redundant. Treat these choices as training-only
hyperparameters, selected before inspecting test outcomes; raw defaults are kept
for compatibility, not as a claim of superior performance.

## Persist the representation

```python
import json
import pickle

features = model.steps[0][1]
feature_names = features.get_feature_names_out()
metadata = features.get_feature_metadata()
with open('spectral_features.json', 'w') as stream:
    json.dump(metadata, stream, indent=2)
with open('pipeline.pkl', 'wb') as stream:
    pickle.dump(model, stream)
```

The metadata records the fit-time package version, exact original embedding,
learned channel means/scales, preprocessing flags, FFT grid, selected bins,
requested spectrum types, and ordered feature names. Store it with exported
features plus source data hashes and your train/test split identifiers. The JSON
is descriptive provenance, not a model loader. Persisting the whole fitted
Pipeline retains the classifier and all preprocessing. Record the Python and
scientific-library versions with the model; cross-version pickle compatibility
is not guaranteed. For pickle/estimator cloning, use ordinary dictionaries for
custom embeddings and named module-level selector functions.

References: scikit-learn's [estimator API](https://scikit-learn.org/stable/developers/develop.html)
and [data-leakage guidance](https://scikit-learn.org/stable/common_pitfalls.html#data-leakage).
