# Residue representations

BLOSUM100 remains the default. Kidera factors provide ten residue-level
physicochemical channels, `KF1`–`KF10`, in published order, covering the 20
canonical amino acids. No ambiguous-residue imputation is performed.

```python
from protspect import SequenceSet
from protspect.embeddings import kidera_factors

features = SequenceSet(embedding=kidera_factors(), compute_cross_spectra=False)
features.from_sequence_list(["ACDE", "ACDF"], spectrum_length=8)
table = features.export_datatable()
```

The same mapping works with `ComplexSet` and the optional ML transformer:

```python
from protspect.sklearn import SpectralTransformer

transformer = SpectralTransformer(embedding=kidera_factors(), spectrum_length=8)
X = transformer.fit_transform(["ACDE", "ACDF"])
```

`kidera_factors()` returns fresh ordinary dictionaries for customization and
scikit-learn cloning/pickling. `protspect.embeddings.kidera` exposes the same
numbers as a read-only nested mapping for direct core use. Use the factory
when passing an explicit embedding to scikit-learn estimators.

These are per-residue factor signals, not ten whole-sequence averages. The
existing Fourier transform retains sequence ordering. Published numbers are
used directly at two-decimal precision; normalization, centering, and fitted
standardization remain explicit options. Cross spectra work as for any other
embedding. Unknown residues raise the existing validation error.

## Reference and provenance

Kidera, A., Konishi, Y., Oka, M., Ooi, T., & Scheraga, H. A. (1985).
Statistical analysis of the physical properties of the 20 naturally occurring
amino acids. *Journal of Protein Chemistry*, 4(1), 23–55.
[doi:10.1007/BF01025492](https://doi.org/10.1007/BF01025492).

The ten factors summarize multiple correlated physical properties; interpret
them as statistical factors rather than ten pure experimental measurements.
Their signs and order are retained exactly.

All 200 numerical values were cross-checked between these pinned reference
implementations, without copying their program logic:

- [CRAN Peptides AAdata](https://github.com/cran/Peptides/blob/ee47d5db6575983bafa9173d85e408c2bca9190d/data/AAdata.RData),
  `AAdata$kideraFactors`, KF1–KF10.
- [peptides.py Kidera tables](https://github.com/althonos/peptides.py/tree/eedec6802b9d17d8895fe5686ae9e371f6be551b/peptides/tables/kidera),
  `KF1.csv`–`KF10.csv`.

The original article verifies the citation and method description; its full
numerical table was checked against the implementations above, rather than a
copy of the paywalled article. Tests pin a digest of the CRAN factor table,
selected vectors, and an independently computed direct DFT.

## AAindex1

Load additional scalar residue properties from an official local AAindex1
file. See [AAindex1 usage and provenance](AAINDEX.md).
