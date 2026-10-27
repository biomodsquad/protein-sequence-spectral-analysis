"""Residue-property channels keyed by channel and canonical residue.

Each outer key names a residue-property channel; each inner mapping gives
its value for a residue. Ambiguous and noncanonical residues are not included.
The historical P-to-A transcription error has been corrected from 2 to -2.
All 400 canonical scores match the NCBI reference matrix:
https://ftp.ncbi.nlm.nih.gov/blast/matrices/BLOSUM100
"""

from types import MappingProxyType

b100 = MappingProxyType(
            {
                'A': MappingProxyType(
                    {'A': 8, 'R': -3, 'N': -4, 'D': -5, 'C': -2, 'Q': -2, 'E': -3, 'G': -1, 'H': -4, 'I': -4,
                     'L': -4, 'K': -2, 'M': -3, 'F': -5, 'P': -2, 'S': 1, 'T': -1, 'W': -6, 'Y': -5, 'V': -2}),

                'R': MappingProxyType(
                    {'A': -3, 'R': 10, 'N': -2, 'D': -5, 'C': -8, 'Q': 0, 'E': -2, 'G': -6, 'H': -1, 'I': -7,
                     'L': -6, 'K': 3, 'M': -4, 'F': -6, 'P': -5, 'S': -3, 'T': -3, 'W': -7, 'Y': -5, 'V': -6}),

                'N': MappingProxyType(
                    {'A': -4, 'R': -2, 'N': 11, 'D': 1, 'C': -5, 'Q': -1, 'E': -2, 'G': -2, 'H': 0, 'I': -7,
                     'L': -7, 'K': -1, 'M': -5, 'F': -7, 'P': -5, 'S': 0, 'T': -1, 'W': -8, 'Y': -5, 'V': -7}),

                'D': MappingProxyType(
                    {'A': -5, 'R': -5, 'N': 1, 'D': 10, 'C': -8, 'Q': -2, 'E': 2, 'G': -4, 'H': -3, 'I': -8,
                     'L': -8, 'K': -3, 'M': -8, 'F': -8, 'P': -5, 'S': -2, 'T': -4, 'W': -10, 'Y': -7, 'V': -8}),

                'C': MappingProxyType(
                    {'A': -2, 'R': -8, 'N': -5, 'D': -8, 'C': 14, 'Q': -7, 'E': -9, 'G': -7, 'H': -8, 'I': -3,
                     'L': -5, 'K': -8, 'M': -4, 'F': -4, 'P': -8, 'S': -3, 'T': -3, 'W': -7, 'Y': -6, 'V': -3}),

                'Q': MappingProxyType(
                    {'A': -2, 'R': 0, 'N': -1, 'D': -2, 'C': -7, 'Q': 11, 'E': 2, 'G': -5, 'H': 1, 'I': -6,
                     'L': -5, 'K': 2, 'M': -2, 'F': -6, 'P': -4, 'S': -2, 'T': -3, 'W': -5, 'Y': -4, 'V': -5}),

                'E': MappingProxyType(
                    {'A': -3, 'R': -2, 'N': -2, 'D': 2, 'C': -9, 'Q': 2, 'E': 10, 'G': -6, 'H': -2, 'I': -7,
                     'L': -7, 'K': 0, 'M': -5, 'F': -8, 'P': -4, 'S': -2, 'T': -3, 'W': -8, 'Y': -7, 'V': -5}),

                'G': MappingProxyType(
                    {'A': -1, 'R': -6, 'N': -2, 'D': -4, 'C': -7, 'Q': -5, 'E': -6, 'G': 9, 'H': -6, 'I': -9,
                     'L': -8, 'K': -5, 'M': -7, 'F': -8, 'P': -6, 'S': -2, 'T': -5, 'W': -7, 'Y': -8, 'V': -8}),

                'H': MappingProxyType(
                    {'A': -4, 'R': -1, 'N': 0, 'D': -3, 'C': -8, 'Q': 1, 'E': -2, 'G': -6, 'H': 13, 'I': -7,
                     'L': -6, 'K': -3, 'M': -5, 'F': -4, 'P': -5, 'S': -3, 'T': -4, 'W': -5, 'Y': 1, 'V': -7}),

                'I': MappingProxyType(
                    {'A': -4, 'R': -7, 'N': -7, 'D': -8, 'C': -3, 'Q': -6, 'E': -7, 'G': -9, 'H': -7, 'I': 8,
                     'L': 2, 'K': -6, 'M': 1, 'F': -2, 'P': -7, 'S': -5, 'T': -3, 'W': -6, 'Y': -4, 'V': 4}),

                'L': MappingProxyType(
                    {'A': -4, 'R': -6, 'N': -7, 'D': -8, 'C': -5, 'Q': -5, 'E': -7, 'G': -8, 'H': -6, 'I': 2,
                     'L': 8, 'K': -6, 'M': 3, 'F': 0, 'P': -7, 'S': -6, 'T': -4, 'W': -5, 'Y': -4, 'V': 0}),

                'K': MappingProxyType(
                    {'A': -2, 'R': 3, 'N': -1, 'D': -3, 'C': -8, 'Q': 2, 'E': 0, 'G': -5, 'H': -3, 'I': -6,
                     'L': -6, 'K': 10, 'M': -4, 'F': -6, 'P': -3, 'S': -2, 'T': -3, 'W': -8, 'Y': -5, 'V': -5}),

                'M': MappingProxyType(
                    {'A': -3, 'R': -4, 'N': -5, 'D': -8, 'C': -4, 'Q': -2, 'E': -5, 'G': -7, 'H': -5, 'I': 1,
                     'L': 3, 'K': -4, 'M': 12, 'F': -1, 'P': -5, 'S': -4, 'T': -2, 'W': -4, 'Y': -5, 'V': 0}),

                'F': MappingProxyType(
                    {'A': -5, 'R': -6, 'N': -7, 'D': -8, 'C': -4, 'Q': -6, 'E': -8, 'G': -8, 'H': -4, 'I': -2,
                     'L': 0, 'K': -6, 'M': -1, 'F': 11, 'P': -7, 'S': -5, 'T': -5, 'W': 0, 'Y': 4, 'V': -3}),

                'P': MappingProxyType(
                    {'A': -2, 'R': -5, 'N': -5, 'D': -5, 'C': -8, 'Q': -4, 'E': -4, 'G': -6, 'H': -5, 'I': -7,
                     'L': -7, 'K': -3, 'M': -5, 'F': -7, 'P': 12, 'S': -3, 'T': -4, 'W': -8, 'Y': -7, 'V': -6}),

                'S': MappingProxyType(
                    {'A': 1, 'R': -3, 'N': 0, 'D': -2, 'C': -3, 'Q': -2, 'E': -2, 'G': -2, 'H': -3, 'I': -5,
                     'L': -6, 'K': -2, 'M': -4, 'F': -5, 'P': -3, 'S': 9, 'T': 2, 'W': -7, 'Y': -5, 'V': -4}),

                'T': MappingProxyType(
                    {'A': -1, 'R': -3, 'N': -1, 'D': -4, 'C': -3, 'Q': -3, 'E': -3, 'G': -5, 'H': -4, 'I': -3,
                     'L': -4, 'K': -3, 'M': -2, 'F': -5, 'P': -4, 'S': 2, 'T': 9, 'W': -7, 'Y': -5, 'V': -1}),

                'W': MappingProxyType(
                    {'A': -6, 'R': -7, 'N': -8, 'D': -10, 'C': -7, 'Q': -5, 'E': -8, 'G': -7, 'H': -5, 'I': -6,
                     'L': -5, 'K': -8, 'M': -4, 'F': 0, 'P': -8, 'S': -7, 'T': -7, 'W': 17, 'Y': 2, 'V': -5}),

                'Y': MappingProxyType(
                    {'A': -5, 'R': -5, 'N': -5, 'D': -7, 'C': -6, 'Q': -4, 'E': -7, 'G': -8, 'H': 1, 'I': -4,
                     'L': -4, 'K': -5, 'M': -5, 'F': 4, 'P': -7, 'S': -5, 'T': -5, 'W': 2, 'Y': 12, 'V': -5}),

                'V': MappingProxyType(
                    {'A': -2, 'R': -6, 'N': -7, 'D': -8, 'C': -3, 'Q': -5, 'E': -5, 'G': -8, 'H': -7, 'I': 4,
                     'L': 0, 'K': -5, 'M': 0, 'F': -3, 'P': -6, 'S': -4, 'T': -1, 'W': -5, 'Y': -5, 'V': 8}),
            })


# Kidera et al. (1985), doi:10.1007/BF01025492; see docs/REPRESENTATIONS.md.
# Published factor order KF1–KF10, residue order ARNDCQEGHILKMFPSTWYV.
_KIDERA_ROWS = (
    ('A', (-1.56, -1.67, -0.97, -0.27, -0.93, -0.78, -0.20, -0.08, 0.21, -0.48)),
    ('R', (0.22, 1.27, 1.37, 1.87, -1.70, 0.46, 0.92, -0.39, 0.23, 0.93)),
    ('N', (1.14, -0.07, -0.12, 0.81, 0.18, 0.37, -0.09, 1.23, 1.10, -1.73)),
    ('D', (0.58, -0.22, -1.58, 0.81, -0.92, 0.15, -1.52, 0.47, 0.76, 0.70)),
    ('C', (0.12, -0.89, 0.45, -1.05, -0.71, 2.41, 1.52, -0.69, 1.13, 1.10)),
    ('Q', (-0.47, 0.24, 0.07, 1.10, 1.10, 0.59, 0.84, -0.71, -0.03, -2.33)),
    ('E', (-1.45, 0.19, -1.61, 1.17, -1.31, 0.40, 0.04, 0.38, -0.35, -0.12)),
    ('G', (1.46, -1.96, -0.23, -0.16, 0.10, -0.11, 1.32, 2.36, -1.66, 0.46)),
    ('H', (-0.41, 0.52, -0.28, 0.28, 1.61, 1.01, -1.85, 0.47, 1.13, 1.63)),
    ('I', (-0.73, -0.16, 1.79, -0.77, -0.54, 0.03, -0.83, 0.51, 0.66, -1.78)),
    ('L', (-1.04, 0.00, -0.24, -1.10, -0.55, -2.05, 0.96, -0.76, 0.45, 0.93)),
    ('K', (-0.34, 0.82, -0.23, 1.70, 1.54, -1.62, 1.15, -0.08, -0.48, 0.60)),
    ('M', (-1.40, 0.18, -0.42, -0.73, 2.00, 1.52, 0.26, 0.11, -1.27, 0.27)),
    ('F', (-0.21, 0.98, -0.36, -1.43, 0.22, -0.81, 0.67, 1.10, 1.71, -0.44)),
    ('P', (2.06, -0.33, -1.15, -0.75, 0.88, -0.45, 0.30, -2.30, 0.74, -0.28)),
    ('S', (0.81, -1.08, 0.16, 0.42, -0.21, -0.43, -1.89, -1.15, -0.97, -0.23)),
    ('T', (0.26, -0.70, 1.21, 0.63, -0.10, 0.21, 0.24, -1.15, -0.56, 0.19)),
    ('W', (0.30, 2.10, -0.72, -1.57, -1.16, 0.57, -0.48, -0.40, -2.30, -0.60)),
    ('Y', (1.38, 1.48, 0.80, -0.56, 0.00, -0.68, -0.31, 1.03, -0.05, 0.53)),
    ('V', (-0.74, -0.71, 2.04, -0.40, 0.50, -0.81, -1.07, 0.06, -0.46, 0.65)),
)

kidera = MappingProxyType({
    f"KF{index + 1}": MappingProxyType({aa: values[index] for aa, values in _KIDERA_ROWS})
    for index in range(10)
})


def kidera_factors():
    """Return independent residue mappings for the ten published Kidera factors.

    Channels KF1–KF10 cover the 20 canonical amino acids in published order
    at two-decimal source precision. No extra scaling, averaging, or unknown
    residue imputation is applied. Fresh dictionaries support customization,
    estimator cloning, and serialization without changing the built-in table.

    Reference: Kidera et al. (1985), Journal of Protein Chemistry 4, 23–55,
    doi:10.1007/BF01025492. See docs/REPRESENTATIONS.md for table provenance.
    """
    return {name: dict(values) for name, values in kidera.items()}
