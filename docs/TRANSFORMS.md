# Sequence transforms

These APIs are available in Protspect 0.1.0a2.
Install with `python -m pip install 'Protspect[ml,wavelets]==0.1.0a2'` for every backend.
PyWavelets is optional and imported only when wavelets are requested.

```python
from protspect import transform_sequence
from protspect.embeddings import kidera_factors

result = transform_sequence(
    'GIGKFLHSAKKFGKAFVGEIMNS', kidera_factors(),
    method='wavelet', length=64, wavelet='db2', level=3,
    center=True, normalization='energy',
)
values, names = result.features('real')
print(result.coordinates, result.metadata)
```

The same embedding mappings work with all transforms, including AAindex1
channels and AAindex2/3 row/column adapters. `Sequence`, `SequenceSet`, `ComplexSet`,
and `SequenceTransformer` use one shared engine for all four transforms. Legacy
`SequenceFFT` and `SpectralTransformer` preserve existing Fourier behavior through
compatibility adapters. Primary exports return no repeated closure bin and apply
no implicit selector. See [the unified API and migration guide](UNIFIED_API.md).

| Backend | Definition / coordinates | Main controls |
| --- | --- | --- |
| `fourier` | Shifted DFT; cycles per residue | `length` |
| `z` | Sum x[n] z^(-n), n=0…N−1; z=r exp(2πif) | positive `radius` |
| `mellin` | Sum x[n] (n+1)^(c−1−2πiν), n=0…N−1 | `mellin_c`, `length` |
| `wavelet` | Multilevel discrete wavelet decomposition; ordered A_L, D_L,…,D_1 | `wavelet`, `level`, `mode` |

The Z transform is evaluated on a circle. Radius 1 equals the DFT; other radii
are computed efficiently by FFT of x[n] r^(-n). Radius below 1 exponentially
upweights later residues; radius above 1 downweights them. This changes sensitivity
to sequence direction and can create large values. Nonfinite results are rejected.
This first implementation does not expose arbitrary spiral CZT contours.
See [SciPy's CZT documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.CZT.html)
for the more general construction.

Mellin uses a **finite discrete sum with unit-spacing quadrature weights**,
not logarithmic resampling or an exact continuous transform. Positive residue
positions start at 1 so log(0) never occurs. The common ν grid is
`fftshift(fftfreq(length, d=log(length)/(length-1)))`; for length 1 use d=1.
This is a selected sampling grid, not a Nyquist claim for nonuniform log positions.
`c=0` supplies reciprocal-position weighting. Length fixes output sampling, but
no fictional padded residues contribute to the sum. Complexity is O(N*length)
per channel. Origin, contour, and discretization matter: there is no general
scale/length invariance guarantee for finite peptide sequences. A future
log-resampled fast Mellin backend should be separately named and benchmarked.

All backends center and normalize over actual residues **before** padding.
Fourier/Z/wavelet then right-pad to `length`. Wavelet boundary extension is applied
by PyWavelets to that padded signal, so both padding and extension affect edge
coefficients. Default `db2`, `mode='symmetric'`; level defaults to
`dwt_max_level(length, filter_length)` and excessive levels are rejected.
Short grids may permit only level 0 (no decomposition). Wavelet output width
can exceed `length`, depending on extension and filter length. Consult
[PyWavelets' DWT conventions](https://pywavelets.readthedocs.io/en/stable/ref/dwt-discrete-wavelet-transform.html).

Use `features('real')` for signed wavelet coefficients; amplitude discards signs.
Phase is defined mathematically but offers little information for real wavelets.
Coefficient indices are accompanied by band and within-band coordinates. They
must not be interpreted as uniformly spaced Fourier frequencies. All backends
support aligned cross-channel and paired-chain coefficient products through
`features(kinds=('cross_coefficients',))`, `transform_pair`, and `ComplexSet`.
Coordinate-dictionary selectors keep band/position and log-frequency meanings
explicit; historical frequency selectors remain in Fourier compatibility exports.

## Exploration notebook

Open [alternative_transforms.ipynb](../examples/alternative_transforms.ipynb)
to edit transform settings, inspect coefficient plots, and compare sequence
direction with the same channels. Run from this checkout using a Python kernel
with `Protspect[ml,wavelets]`, matplotlib, pandas, and nbformat installed.

## Fixed-width machine-learning features

```python
from protspect.sklearn import SequenceTransformer
from sklearn.pipeline import make_pipeline
from sklearn.linear_model import LogisticRegression

features = SequenceTransformer(
    method='wavelet', length=64, component='real',
    center=True, normalization='energy',
    transform_options={'wavelet': 'haar', 'level': 4, 'mode': 'periodization'},
)
model = make_pipeline(features, LogisticRegression())
```

Use `method='fourier'`, `'z'`, or `'mellin'` with the same estimator.
`transform_options={'radius': 1.05}` configures Z;
`transform_options={'mellin_c': 0.0}` configures Mellin. Grid-search the options
mapping as a whole (`transform_options=[{...}, {...}]`). Fit inside each CV fold.
When length is omitted, fit chooses the next power of two covering training
sequences only. Longer held-out sequences raise; they never expand the schema.
Set an explicit grid if length coverage is known in advance.
Optional `standardize_channels=True` learns training-residue moments only.
Snapshots protect fitted behavior from later constructor-dictionary mutations;
call fit again to apply changes. Features contain all coefficients by default,
including redundant Fourier conjugates. Compare representations on held-out
peptides before claiming improvement; representation-atlas proximity alone does
not measure transform performance.
