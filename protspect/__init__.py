"""Protspect: spectral features for protein and peptide sequences."""

from .features import ComplexSet, SequenceFFT, SequenceSet

from .transforms import TransformResult, transform_sequence

__all__ = ["ComplexSet", "SequenceFFT", "SequenceSet", "TransformResult", "transform_sequence"]
__version__ = "0.1.0a1"
