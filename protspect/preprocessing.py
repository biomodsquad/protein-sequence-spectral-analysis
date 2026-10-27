"""Explicit per-sequence preprocessing before padding and transformation."""

import numpy as np

NORMALIZATIONS = ("none", "length", "sqrt_length", "energy")


def validate_preprocessing(normalization, center):
    """Validate normalization and centering without changing caller parameters."""
    if normalization not in NORMALIZATIONS:
        raise ValueError(f"normalization must be one of {NORMALIZATIONS}")
    if not isinstance(center, (bool, np.bool_)):
        raise ValueError("center must be a boolean")


def prepare_signal(vector, normalization, center):
    """Center real residues, then scale, leaving padding out of both operations.

    ``length`` and ``sqrt_length`` divide by the original residue count and
    its square root. ``energy`` divides by the channel's L2 norm *after*
    optional centering, producing unit energy for nonzero signals. Zero-energy
    channels remain zero. None of these options is the FFT's orthogonal norm.
    """
    values = np.asarray(vector, dtype=float)
    if center:
        values = (np.zeros_like(values) if np.all(values == values[0])
                  else values - np.sum(values / values.size))
    if normalization == "length":
        values = values / values.size
    elif normalization == "sqrt_length":
        values = values / np.sqrt(values.size)
    elif normalization == "energy":
        # Scaling first avoids overflow/underflow when computing the L2 norm.
        peak = np.max(np.abs(values))
        if peak:
            scaled = values / peak
            values = scaled / np.sqrt(np.sum(scaled * scaled))
    if not np.all(np.isfinite(values)):
        raise ValueError("preprocessing produced nonfinite channel values")
    return values
