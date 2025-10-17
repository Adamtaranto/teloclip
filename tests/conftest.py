"""
Pytest configuration and fixtures for teloclip tests.
"""

from pathlib import Path
import tempfile

import pytest


@pytest.fixture
def temp_dir():
    """Create a temporary directory for test files."""
    with tempfile.TemporaryDirectory() as temp_dir:
        yield Path(temp_dir)


@pytest.fixture
def sample_fasta_content():
    """Sample FASTA content for testing."""
    return """>contig01
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>contig02
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG
>contig03
TTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAA
TTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAA
TTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAA
TTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAA
"""


@pytest.fixture
def sample_fai_content():
    """Sample FAI index content for testing."""
    return """contig01\t120\t9\t60\t61
contig02\t180\t139\t60\t61
contig03\t240\t329\t60\t61
"""


@pytest.fixture
def sample_sam_header():
    """Sample SAM header for testing."""
    return """@HD\tVN:1.3\tSO:coordinate
@SQ\tSN:contig01\tLN:120
@SQ\tSN:contig02\tLN:180
@SQ\tSN:contig03\tLN:240
@PG\tID:minimap2\tPN:minimap2\tVN:2.10-r761\tCL:minimap2 -x sr -a test.fna reads.fq.gz
"""


@pytest.fixture
def sample_sam_alignments():
    """Sample SAM alignments with various clip patterns."""
    return [
        # Left terminal clip - alignment starts at position 1 with left soft clip
        'read01\t0\tcontig01\t1\t60\t20S100M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t*',
        # Right terminal clip - alignment ends near contig end with right soft clip
        'read02\t0\tcontig01\t21\t60\t80M20S\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t*',
        # Both ends clipped
        'read03\t0\tcontig02\t1\t60\t15S150M15S\t*\t0\t0\tGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG\t*',
        # No clips - should not be output
        'read04\t0\tcontig01\t10\t60\t100M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t*',
        # Complex CIGAR with insertions/deletions and terminal clip
        'read05\t0\tcontig03\t201\t60\t30M5I10M3D20M25S\t*\t0\t0\tTTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTA\t*',
        # Low anchor - should be filtered with high min_anchor
        'read06\t0\tcontig01\t110\t60\t10M30S\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t*',
    ]


@pytest.fixture
def telomeric_sequences():
    """Sample sequences with telomeric repeats for motif testing."""
    return {
        'telomeric': 'TTAGGGTTAGGGTTAGGGTTAGGGTTAGGG',
        'reverse_telomeric': 'CCCTAACCCTAACCCTAACCCTAACCCTAA',
        'mixed': 'ATCGTTAGGGCCCTAATCGTTAGGGCCCTAA',
        'no_telomeric': 'ATCGATCGATCGATCGATCGATCGATCGATCG',
        'fuzzy_telomeric': 'TTAGGTTTAGGGTTTAGGGTTTAGGG',  # For fuzzy matching tests
    }


@pytest.fixture
def create_test_files(temp_dir):
    """Create test files in temporary directory."""

    def _create_files(fasta_content=None, fai_content=None, sam_content=None):
        files = {}

        if fasta_content:
            fasta_file = temp_dir / 'test.fasta'
            fasta_file.write_text(fasta_content)
            files['fasta'] = str(fasta_file)

        if fai_content:
            fai_file = temp_dir / 'test.fasta.fai'
            fai_file.write_text(fai_content)
            files['fai'] = str(fai_file)

        if sam_content:
            sam_file = temp_dir / 'test.sam'
            sam_file.write_text(sam_content)
            files['sam'] = str(sam_file)

        return files

    return _create_files


@pytest.fixture
def cigar_test_cases():
    """Test cases for CIGAR string parsing."""
    return [
        # (CIGAR, expected_left_clip, expected_right_clip, expected_aligned_bases, expected_ref_length)
        ('100M', None, None, 100, 100),
        ('20S80M', 20, None, 80, 80),
        ('80M20S', None, 20, 80, 80),
        ('10S60M30S', 10, 30, 60, 60),
        ('50M10I40M', None, None, 90, 90),
        ('50M5D40M', None, None, 90, 95),
        ('30M10I5D20M15S', None, 15, 50, 55),
        ('15H20S30M40=20X10S5H', 20, 10, 90, 90),
        ('100S', 100, None, 0, 0),
        ('25M', None, None, 25, 25),
    ]
