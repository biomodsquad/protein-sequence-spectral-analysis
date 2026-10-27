"""Offline AAindex records and conversion to Protspect embedding channels.

No database is bundled or fetched during import, parsing, or feature extraction.
Use download_aaindex1/2/3 explicitly to retrieve official files, then from_file
for offline use. Matrix conversion requires an explicit profile axis.
"""

from collections.abc import Mapping
from dataclasses import dataclass
import hashlib
import math
import os
from pathlib import Path
import re
import tempfile
from types import MappingProxyType
from urllib.request import urlopen

AAINDEX1_URL = 'https://www.genome.jp/ftp/db/community/aaindex/aaindex1'
AAINDEX2_URL = 'https://www.genome.jp/ftp/db/community/aaindex/aaindex2'
AAINDEX3_URL = 'https://www.genome.jp/ftp/db/community/aaindex/aaindex3'
_ALPHABET = 'ARNDCQEGHILKMFPSTWYV'
_ACCESSION = re.compile(r'[A-Z]{4}[0-9]{6}')


def _check_digest(content, expected):
    """Return the raw-file digest, rejecting malformed or mismatched pins."""
    digest = hashlib.sha256(content).hexdigest()
    if expected is not None:
        if not isinstance(expected, str) or not re.fullmatch(r'[0-9a-fA-F]{64}', expected):
            raise ValueError('expected_sha256 must be a 64-character hexadecimal digest')
        if digest != expected.lower():
            raise ValueError(f'AAindex SHA256 mismatch: expected {expected.lower()}, got {digest}')
    return digest


@dataclass(frozen=True)
class AAIndexEntry:
    """One AAindex1 property and its literature metadata.

    scores follows ARNDCQEGHILKMFPSTWYV, with None for published NA values.
    values exposes a read-only residue mapping; metadata preserves source text
    without interpreting units or treating correlations as interchangeable indices.
    """

    accession: str
    description: str
    scores: tuple
    authors: str = ''
    title: str = ''
    journal: str = ''
    references: str = ''
    comments: str = ''
    correlations: str = ''

    @property
    def values(self):
        """Read-only canonical residue values, including explicit missing values."""
        return MappingProxyType(dict(zip(_ALPHABET, self.scores)))

    @property
    def missing_residues(self):
        """Canonical residues lacking a published value, in canonical order."""
        return tuple(aa for aa, value in zip(_ALPHABET, self.scores) if value is None)

    @property
    def complete(self):
        """Whether all twenty canonical residues have finite published values."""
        return not self.missing_residues


def _record(fields):
    """Decode declared two-row residue order instead of guessing alphabetical order."""
    accession = ' '.join(fields.get('H', []))
    if not _ACCESSION.fullmatch(accession):
        raise ValueError(f'invalid AAindex1 accession: {accession!r}')
    if 'M' in fields:
        raise ValueError(f'{accession}: matrix records require AAIndex2 or AAIndex3')
    if 'I' not in fields or not fields.get('D'):
        raise ValueError(f'{accession}: AAindex1 requires description and I data fields')
    lines = fields['I']
    pairs = [pair.split('/') for pair in lines[0].split()]
    if len(pairs) != 10 or any(len(pair) != 2 for pair in pairs):
        raise ValueError(f'{accession}: invalid I residue-order header')
    order = [pair[0] for pair in pairs] + [pair[1] for pair in pairs]
    if len(set(order)) != 20 or set(order) != set(_ALPHABET):
        raise ValueError(f'{accession}: I header must contain exactly the 20 canonical residues')
    tokens = ' '.join(lines[1:]).split()
    if len(tokens) != 20:
        raise ValueError(f'{accession}: expected 20 index values, got {len(tokens)}')
    values = {}
    for residue, token in zip(order, tokens):
        try:
            value = None if token == 'NA' else float(token)
        except ValueError as error:
            raise ValueError(f'{accession}: invalid value for {residue}: {token!r}') from error
        if value is not None and not math.isfinite(value):
            raise ValueError(f'{accession}: nonfinite value for {residue}')
        values[residue] = value
    def text(key):
        return ' '.join(fields.get(key, []))
    return AAIndexEntry(accession, text('D'), tuple(values[aa] for aa in _ALPHABET),
                        text('A'), text('T'), text('J'), text('R'), text('*'), text('C'))


def _parse(text, record_parser=_record):
    """Parse tagged records, continuation lines, and exact // terminators."""
    entries, fields, active = {}, {}, None
    for number, line in enumerate(text.splitlines(), 1):
        if not line.strip():
            continue
        if line.strip() == '//':
            if not fields:
                raise ValueError(f'line {number}: empty AAindex record')
            entry = record_parser(fields)
            if entry.accession in entries:
                raise ValueError(f'duplicate accession: {entry.accession}')
            entries[entry.accession] = entry
            fields, active = {}, None
        elif line[0].isspace():
            if active is None:
                raise ValueError(f'line {number}: continuation without a field')
            fields[active].append(line.strip())
        elif line[0] in 'HDRATJ*CIM' and (len(line) == 1 or line[1].isspace()):
            active = line[0]
            if active in fields and active != '*':
                raise ValueError(f'line {number}: duplicate {active} field')
            fields.setdefault(active, []).append(line[2:].strip())
        else:
            raise ValueError(f'line {number}: unrecognized AAindex field')
    if fields:
        raise ValueError('unterminated AAindex record: expected //')
    if not entries:
        raise ValueError('AAindex file contains no records')
    return entries


class AAIndex(Mapping):
    """Read-only catalogue of AAindex1 records; load files explicitly with from_file.

    Entries may contain missing values; embedding conversion requires complete
    entries. No implicit scaling, imputation, filtering, or network access occurs.
    """

    _record_parser = staticmethod(_record)

    def __init__(self, entries, *, sha256, source_path, source_url=None, version=None):
        self._entries = MappingProxyType(dict(entries))
        self.sha256 = sha256
        self.source_path = source_path
        self.source_url = source_url
        self.version = version

    @classmethod
    def from_file(cls, path, *, expected_sha256=None, source_url=None, version=None):
        """Read a local UTF-8 file for this catalogue section, optionally pinning SHA256.

        source_url and version are caller-supplied provenance labels: AAindex
        records do not declare a database version. The source path and actual
        raw-byte checksum are always recorded. This method never downloads data.
        """
        path = Path(path).resolve()
        content = path.read_bytes()
        digest = _check_digest(content, expected_sha256)
        entries = _parse(content.decode('utf-8-sig'), cls._record_parser)
        return cls(entries, sha256=digest, source_path=str(path),
                   source_url=source_url, version=version)

    def __getitem__(self, accession):
        """Return a record by exact accession, raising KeyError if absent."""
        return self._entries[accession]

    def __iter__(self):
        """Iterate accessions in file order."""
        return iter(self._entries)

    def __len__(self):
        """Number of parsed records, including incomplete entries."""
        return len(self._entries)

    def search(self, query, *, complete_only=False):
        """Return records matching a case-insensitive substring of citation or description."""
        if not isinstance(query, str) or not query.strip():
            raise ValueError('query must be a nonempty string')
        query = query.strip().casefold()
        return tuple(entry for entry in self.values()
                     if (not complete_only or entry.complete) and query in ' '.join(
                         (entry.accession, entry.description, entry.authors, entry.title,
                          entry.journal, entry.references, entry.comments)).casefold())

    def _select(self, accessions):
        """Keep explicit caller order and reject ambiguous duplicate channel requests."""
        ids = (accessions,) if isinstance(accessions, str) else tuple(accessions)
        if not ids or any(not isinstance(key, str) for key in ids):
            raise ValueError('select at least one accession string')
        if len(set(ids)) != len(ids):
            raise ValueError('duplicate accessions would produce duplicate channels')
        return [self[key] for key in ids]

    def embedding(self, accessions):
        """Return fresh channel-to-residue dictionaries in requested accession order.

        Each accession names one channel. Published values are unchanged.
        Incomplete entries raise ValueError listing the missing residues; no
        imputation is performed. Unknown accessions raise KeyError. Returned
        mappings support Protspect core classes and scikit-learn estimators.
        """
        entries = self._select(accessions)
        missing = [f'{e.accession}: {", ".join(e.missing_residues)}'
                   for e in entries if not e.complete]
        if missing:
            raise ValueError('AAindex1 entries have missing values: ' + '; '.join(missing))
        return {entry.accession: dict(entry.values) for entry in entries}

    def metadata(self, accessions):
        """Return JSON-compatible provenance and full selected source values.

        Export this together with fitted preprocessing metadata. Missing source
        values remain None (JSON null); incomplete entries can be inspected
        without converting them into an embedding.
        """
        entries = self._select(accessions)
        return {'database': 'AAindex1', 'version': self.version,
                'source_path': self.source_path, 'source_url': self.source_url,
                'sha256': self.sha256, 'residue_order': _ALPHABET,
                'channels': [e.accession for e in entries],
                'entries': [dict(accession=e.accession, description=e.description,
                                 authors=e.authors, title=e.title, journal=e.journal,
                                 references=e.references, comments=e.comments,
                                 correlations=e.correlations, values=dict(e.values))
                            for e in entries]}


def _download_database(path, *, expected_sha256=None, timeout=30,
                       catalogue_class=AAIndex, url=AAINDEX1_URL):
    """Download and validate a selected official section into a new cache file.

    Validate checksum and all records before writing. Refuse to overwrite any
    existing file; reuse it offline with its catalogue from_file. Returns the parsed
    catalogue with the official source URL. No third-party HTTP dependency is
    required. Downloaded data is separate from Protspect's package license.
    """
    path = Path(path).resolve()
    if path.exists():
        raise FileExistsError(
            f'{path} exists; load it with {catalogue_class.__name__}.from_file')
    if not math.isfinite(timeout) or timeout <= 0:
        raise ValueError('timeout must be a positive finite number')
    if expected_sha256 is not None and (not isinstance(expected_sha256, str)
            or not re.fullmatch(r'[0-9a-fA-F]{64}', expected_sha256)):
        raise ValueError('expected_sha256 must be a 64-character hexadecimal digest')
    with urlopen(url, timeout=timeout) as response:
        content = response.read()
    digest = _check_digest(content, expected_sha256)
    entries = _parse(content.decode('utf-8-sig'), catalogue_class._record_parser)
    path.parent.mkdir(parents=True, exist_ok=True)
    # A same-directory temporary file keeps malformed/partial downloads out of
    # the cache; hard-link creation atomically refuses a concurrently created path.
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as stream:
            temporary = Path(stream.name)
            stream.write(content)
        os.link(temporary, path)
    finally:
        if temporary is not None:
            temporary.unlink()
    return catalogue_class(entries, sha256=digest, source_path=str(path), source_url=url)


@dataclass(frozen=True)
class AAIndexMatrixEntry:
    """An AAindex2/3 matrix with original axes and literature metadata.

    scores is an immutable tuple of row tuples. NA and standalone - remain None. Only
    lower-triangle input is reflected to reconstruct symmetry; full matrices
    retain direction and rectangular shape, including extra residue symbols.
    """

    accession: str
    description: str
    rows: str
    columns: str
    scores: tuple
    storage: str
    authors: str = ''
    title: str = ''
    journal: str = ''
    references: str = ''
    comments: str = ''
    correlations: str = ''

    @property
    def values(self):
        """Read-only mappings indexed first by row residue, then column residue."""
        return MappingProxyType({aa: MappingProxyType(dict(zip(self.columns, row)))
                                 for aa, row in zip(self.rows, self.scores)})

    @property
    def shape(self):
        """Original row and column dimensions, including noncanonical symbols."""
        return len(self.rows), len(self.columns)

    @property
    def missing_pairs(self):
        """Ordered (row, column) pairs with missing source scores."""
        return tuple((r, c) for r, row in zip(self.rows, self.scores)
                     for c, value in zip(self.columns, row) if value is None)

    @property
    def complete(self):
        """Whether every cell of the full native matrix is present."""
        return not self.missing_pairs

    @property
    def canonical_complete(self):
        """Whether the canonical 20-by-20 submatrix exists and is complete."""
        return (set(_ALPHABET) <= set(self.rows) and set(_ALPHABET) <= set(self.columns)
                and not any(r in _ALPHABET and c in _ALPHABET for r, c in self.missing_pairs))


def _matrix_record(fields):
    """Parse full or lower-triangular matrices using the declared axes."""
    def text(key):
        return ' '.join(fields.get(key, []))

    accession = text('H')
    if not re.fullmatch(r'[A-Z_]{4}[0-9]{6}', accession):
        raise ValueError(f'invalid AAindex matrix accession: {accession!r}')
    if 'I' in fields or 'M' not in fields or not fields.get('D'):
        raise ValueError(f'{accession}: matrix catalogue requires M data and description')
    lines = fields['M']
    header = re.fullmatch(r'rows\s*=\s*([A-Za-z*\-]+)\s*,\s*cols\s*=\s*([A-Za-z*\-]+)',
                          lines[0])
    if header is None:
        raise ValueError(f'{accession}: invalid M axis header')
    rows, columns = header.groups()
    if any(len(set(axis)) != len(axis) for axis in (rows, columns)):
        raise ValueError(f'{accession}: matrix axes contain duplicate symbols')
    tokens = ' '.join(lines[1:]).split()
    values = []
    for token in tokens:
        try:
            value = None if token in ('NA', '-') else float(token)
        except ValueError as error:
            raise ValueError(f'{accession}: invalid matrix value: {token!r}') from error
        if value is not None and not math.isfinite(value):
            raise ValueError(f'{accession}: nonfinite matrix value')
        values.append(value)
    nrows, ncols = len(rows), len(columns)
    if len(values) == nrows*ncols:
        scores = tuple(tuple(values[i*ncols:(i+1)*ncols]) for i in range(nrows))
        storage = 'full'
    elif (rows == columns and len(values) == nrows*(nrows+1)//2
          and [len(line.split()) for line in lines[1:]] == list(range(1, nrows+1))):
        matrix = [[None]*ncols for _ in rows]
        cursor = 0
        for i in range(nrows):
            for j in range(i+1):
                matrix[i][j] = matrix[j][i] = values[cursor]
                cursor += 1
        scores = tuple(tuple(row) for row in matrix)
        storage = 'lower_triangle'
    else:
        raise ValueError(f'{accession}: {len(values)} matrix values do not match declared '
                         f'{nrows}x{ncols} axes (full or identical-axis lower triangle)')
    return AAIndexMatrixEntry(accession, text('D'), rows, columns, scores, storage,
                              text('A'), text('T'), text('J'), text('R'), text('*'), text('C'))


class _MatrixCatalogue(AAIndex):
    """Shared matrix catalogue behavior; use AAIndex2 or AAIndex3."""

    _record_parser = staticmethod(_matrix_record)

    def embedding(self, accessions, *, axis):
        """Convert selected canonical submatrices into explicit row/column profiles.

        axis='rows': channel '<ID>_rows_A' maps residue R to M[A,R].
        axis='columns': channel '<ID>_columns_A' maps residue R to M[R,A].
        Every matrix produces twenty channels in ARNDCQEGHILKMFPSTWYV order,
        each covering that same canonical input alphabet. Extra symbols remain
        available in the catalogue but are omitted from this canonical adapter.
        Incomplete canonical blocks raise ValueError. No averaging, scaling,
        imputation, or symmetrization of full matrices is performed. Returned
        ordinary dictionaries support estimator cloning and serialization.
        """
        if axis not in ('rows', 'columns'):
            raise ValueError("axis must be explicitly 'rows' or 'columns'")
        entries = self._select(accessions)
        result = {}
        for entry in entries:
            absent = {'rows': ''.join(aa for aa in _ALPHABET if aa not in entry.rows),
                      'columns': ''.join(aa for aa in _ALPHABET if aa not in entry.columns)}
            if any(absent.values()):
                raise ValueError(f'{entry.accession}: canonical residues absent from axes: {absent}')
            missing = [(r, c) for r, c in entry.missing_pairs if r in _ALPHABET and c in _ALPHABET]
            if missing:
                raise ValueError(f'{entry.accession}: missing canonical matrix values: {missing}')
            matrix = entry.values
            for channel in _ALPHABET:
                result[f'{entry.accession}_{axis}_{channel}'] = {
                    aa: matrix[channel][aa] if axis == 'rows' else matrix[aa][channel]
                    for aa in _ALPHABET}
        return result

    def metadata(self, accessions, *, axis=None):
        """Export source matrices and optionally the exact canonical profile schema.

        With axis=None, incomplete/native records remain inspectable. Supply
        the same explicit axis used for embedding() to record derived channel
        names, orientation, and canonical-submatrix policy alongside the source.
        Full axes and expanded values are always retained with the raw-file hash.
        """
        entries = self._select(accessions)
        channels = None if axis is None else list(self.embedding(
            [entry.accession for entry in entries], axis=axis))
        return {'database': self.database_name, 'version': self.version,
                'source_path': self.source_path, 'source_url': self.source_url,
                'sha256': self.sha256,
                'representation': None if axis is None else {
                    'kind': 'canonical_matrix_profiles', 'axis': axis,
                    'residue_order': _ALPHABET, 'extra_symbols': 'excluded',
                    'normalization': 'none', 'imputation': 'none'},
                'channels': channels,
                'entries': [dict(accession=e.accession, description=e.description,
                                 authors=e.authors, title=e.title, journal=e.journal,
                                 references=e.references, comments=e.comments,
                                 correlations=e.correlations, rows=e.rows, columns=e.columns,
                                 storage=e.storage,
                                 values={r: dict(row) for r, row in e.values.items()})
                            for e in entries]}


class AAIndex2(_MatrixCatalogue):
    """Offline AAindex2 substitution matrix catalogue, loaded with from_file()."""

    database_name = 'AAindex2'


class AAIndex3(_MatrixCatalogue):
    """Offline AAindex3 contact-potential catalogue, loaded with from_file().

    Matrix-derived profile channels are sequence representations. Computing
    their Fourier spectra does not infer structural contacts or contact energies.
    """

    database_name = 'AAindex3'


def download_aaindex1(path, *, expected_sha256=None, timeout=30):
    """Download/validate official AAindex1 into a new file; return an AAIndex catalogue.

    Never overwrite an existing file. Reuse it offline with AAIndex.from_file.
    An optional SHA256 pin enforces exact source bytes before any cache write.
    """
    return _download_database(path, expected_sha256=expected_sha256, timeout=timeout)


def download_aaindex2(path, *, expected_sha256=None, timeout=30):
    """Download/validate official AAindex2 into a new file; return an AAIndex2 catalogue.

    Existing files are never overwritten. Reuse snapshots offline with from_file.
    """
    return _download_database(path, expected_sha256=expected_sha256, timeout=timeout,
                              catalogue_class=AAIndex2, url=AAINDEX2_URL)


def download_aaindex3(path, *, expected_sha256=None, timeout=30):
    """Download/validate official AAindex3 into a new file; return an AAIndex3 catalogue.

    Existing files are never overwritten. Reuse snapshots offline with from_file.
    """
    return _download_database(path, expected_sha256=expected_sha256, timeout=timeout,
                              catalogue_class=AAIndex3, url=AAINDEX3_URL)
