"""
Benchmarks for the canonical overhang predicate.

``teloclip.core.overhang.classify`` runs once per alignment in every sub-command, so
it sits on the hottest path in the codebase. Its cost is dominated by CIGAR
parsing, which is regex-driven, so the parsing helpers are measured separately
to make a regression attributable.
"""

import pytest

from teloclip.core.overhang import (
    anchor_length,
    cigar_ops_from_string,
    cigar_ops_from_tuples,
    classify,
    ends_from_sam_fields,
    reference_span,
)

CONTIG_LENGTH = 1_000_000


@pytest.mark.benchmark
def test_cigar_parsing(cigar_strings):
    """Benchmark CIGAR string parsing across representative shapes."""
    for _ in range(1000):
        for cigar in cigar_strings:
            cigar_ops_from_string(cigar)


@pytest.mark.benchmark
def test_cigar_from_pysam_tuples():
    """Benchmark the pysam CIGAR conversion used by the extend path."""
    tuples = [(4, 50), (0, 100), (2, 5), (0, 95), (1, 3), (4, 20)]
    for _ in range(8000):
        cigar_ops_from_tuples(tuples)


@pytest.mark.benchmark
def test_span_and_anchor(cigar_strings):
    """Benchmark the reference-span and anchor-length reductions."""
    parsed = [cigar_ops_from_string(c) for c in cigar_strings]
    for _ in range(1000):
        for ops in parsed:
            reference_span(ops)
            anchor_length(ops)


@pytest.mark.benchmark
def test_classify_from_sam_fields(sam_lines):
    """Benchmark the full per-alignment path: parse fields, then classify.

    This is what filter and extract do for every record, so it is the closest
    proxy for their throughput.
    """
    records = [
        line.rstrip('\n').split('\t') for line in sam_lines if not line.startswith('@')
    ]

    for fields in records:
        ends = ends_from_sam_fields(fields, CONTIG_LENGTH)
        classify(ends, max_break=50, min_clip=1, min_anchor=100)
