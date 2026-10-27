"""Frequency selectors for spectra measured in cycles per residue."""

import numpy as np


def frequency_select_physical(freq):
    """Select DC and positive dyadic frequencies (1/2, 1/4, ...).

    Parameters
    ----------
    freq : array-like of float
        One-dimensional frequency grid in cycles per residue.

    Returns
    -------
    numpy.ndarray of bool
        Mask selecting grid points equal, within floating-point tolerance, to
        zero or reciprocal powers of two resolvable on the grid. For grids
        without dyadic bins, only DC may be selected.
    """
    freq = np.asarray(freq, dtype=float)
    if freq.ndim != 1 or not np.all(np.isfinite(freq)):
        raise ValueError("freq must be a finite one-dimensional frequency grid")
    if freq.size == 0:
        return np.zeros(0, dtype=bool)
    positive = freq[freq > 0]
    mask = freq == 0
    if positive.size:
        count = max(1, int(np.ceil(-np.log2(positive.min()))))
        for exponent in range(1, count + 1):
            mask |= np.isclose(freq, 2.0 ** -exponent, rtol=1e-10, atol=0)
    return mask


def frequency_select_positive(freq):
    """Return a boolean mask selecting nonnegative frequencies, including DC."""
    return np.asarray(freq) >= 0


def frequency_selelect_positive(freq):
    """Compatibility alias for the historically misspelled positive selector."""
    return frequency_select_positive(freq)
