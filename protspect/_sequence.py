"""Shared input validation and embedding; independent of transform backends."""

from collections.abc import Mapping
from numbers import Integral

import numpy as np


def _sequences(sequences):
    """Materialize and validate a nonempty iterable of nonempty strings."""
    if isinstance(sequences, str):
        raise TypeError("sequences must be an iterable of strings, not a string")
    result = list(sequences)
    if not result or any(not isinstance(s, str) or not s for s in result):
        raise ValueError("provide at least one nonempty sequence string")
    return result


def _length(length, minimum):
    """Choose a power-of-two length or validate an explicit FFT length."""
    if length is None:
        return 1 << (minimum - 1).bit_length()
    if isinstance(length, bool) or not isinstance(length, Integral) or length < minimum:
        raise ValueError("spectrum length must be an integer at least as long as every sequence")
    return int(length)


def _attributes(attributes, count):
    """Convert column-oriented metadata to rows without silently truncating."""
    if attributes is None:
        return [None] * count
    if not isinstance(attributes, Mapping):
        raise TypeError("attributes must be a mapping of names to columns")
    columns = {name: list(values) for name, values in attributes.items()}
    if any(len(values) != count for values in columns.values()):
        raise ValueError("every attribute column must match the number of sequences")
    return [{name: values[i] for name, values in columns.items()} for i in range(count)]


def embed_vectors(sequence, embeddings):
    """Validate a residue mapping and return independent finite real channels."""
    if not isinstance(embeddings, Mapping) or not embeddings:
        raise ValueError("embeddings must be a nonempty mapping of channels")
    vectors = {}
    for name, values in embeddings.items():
        if not isinstance(name, str) or not name:
            raise ValueError("embedding channel names must be nonempty strings")
        if not isinstance(values, Mapping):
            raise ValueError("embedding channels must be residue mappings")
        try:
            scores = [values[aa] for aa in sequence]
            if any(np.ndim(v) != 0 or np.iscomplexobj(v) for v in scores):
                raise ValueError("embedding values must be finite real scalars")
            vector = np.asarray(scores, dtype=float)
        except KeyError as error:
            raise ValueError(f"missing residue {error.args[0]!r} in channel {name!r}") from error
        if vector.ndim != 1 or not np.all(np.isfinite(vector)):
            raise ValueError("embedding values must be finite real scalars")
        vectors[name] = vector
    return vectors
