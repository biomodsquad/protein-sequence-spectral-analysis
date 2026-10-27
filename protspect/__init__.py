"""Protspect: spectral features for protein and peptide sequences."""

from .features import ComplexSet, Sequence, SequenceFFT, SequenceSet

from .transforms import TransformResult, transform_sequence, transform_pair

__all__ = ["Sequence", "ComplexSet", "SequenceFFT", "SequenceSet", "TransformResult", "transform_sequence", "transform_pair"]
__version__ = "0.1.0a1"
