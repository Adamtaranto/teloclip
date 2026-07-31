"""
Shared fixtures for the teloclip benchmark suite.

Benchmarks live outside ``tests/`` so they are not collected by the normal test
run, which is configured with coverage instrumentation that would distort the
measurements.

Inputs are synthesised at a size where the measurement reflects the algorithm
rather than fixture overhead. The mock dataset in ``tests/data`` is only a few
hundred records, which is too small to be informative on its own, so the SAM
benchmarks scale it up.
"""

from pathlib import Path

import pytest

DATA_DIR = Path(__file__).resolve().parents[1] / 'tests' / 'data'
TEST_SAM = DATA_DIR / 'test.sam'
TEST_BAM = DATA_DIR / 'test.bam'
TEST_FNA = DATA_DIR / 'test.fna'
TEST_FAI = DATA_DIR / 'test.fna.fai'

# Contigs in the mock dataset, as name -> length.
CONTIGS = {
    'contig01': 560,
    'contig02': 1610,
    'contig03': 2310,
    'contig04': 906,
}


@pytest.fixture(scope='session')
def sam_lines():
    """
    Provide SAM records repeated to a size worth measuring.

    Returns
    -------
    list of str
        Header lines followed by roughly 25,000 alignment records.
    """
    raw = TEST_SAM.read_text().splitlines(keepends=True)
    headers = [line for line in raw if line.startswith('@')]
    records = [line for line in raw if not line.startswith('@')]

    # Repeat to a stable, meaningful record count regardless of the fixture size.
    repeats = max(1, 25_000 // max(1, len(records)))
    return headers + records * repeats


@pytest.fixture(scope='session')
def cigar_strings():
    """
    Provide a spread of CIGAR strings covering the shapes teloclip encounters.

    Returns
    -------
    list of str
        Clipped, unclipped, indel-bearing and both-ended CIGARs.
    """
    return [
        '50S100M',
        '100M50S',
        '20S200M30S',
        '150M',
        '30M5D10M5I105M20S',
        '10S30M5N60M2X8=15S',
        '250M',
        '75S175M',
    ]


@pytest.fixture(scope='session')
def contigs():
    """
    Provide the mock assembly's contig lengths.

    Returns
    -------
    dict
        Mapping of contig name to length.
    """
    return dict(CONTIGS)


@pytest.fixture(scope='session')
def bam_path():
    """
    Provide the path to the indexed mock BAM.

    Returns
    -------
    Path
        Location of ``tests/data/test.bam``.
    """
    return TEST_BAM


@pytest.fixture(scope='session')
def fasta_path():
    """
    Provide the path to the mock reference FASTA.

    Returns
    -------
    Path
        Location of ``tests/data/test.fna``.
    """
    return TEST_FNA
