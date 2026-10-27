"""Small, self-contained example requiring only the package dependencies."""

from protspect import ComplexSet, SequenceSet

sequences = SequenceSet(compute_cross_spectra=False)
sequences.from_sequence_list(
    ["ACDE", "ACDF"],
    spectrum_length=8,
    attributes={"name": ["peptide-1", "peptide-2"]},
)
print(sequences.export_datatable(attributes=["name"]))

pairs = ComplexSet()
pairs.from_sequence_list(["ACDE"], ["ACDF"], spectrum_length=8)
print(pairs.export_datatable().shape)
