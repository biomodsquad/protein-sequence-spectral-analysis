"""Core-only unified API example; both collection types support every backend."""

from protspect import ComplexSet, SequenceSet

sequences = SequenceSet(method='fourier', compute_cross_spectra=False)
sequences.from_sequence_list(
    ["ACDE", "ACDF"], length=8,
    attributes={"name": ["peptide-1", "peptide-2"]},
)
print(sequences.export_features(attributes=["name"]))

pairs = ComplexSet(method='z', transform_options={'radius': 1.05})
pairs.from_sequence_list(["ACDE"], ["ACDF"], length=8)
print(pairs.export_features(('amplitude', 'phase')).shape)
