"""
Unit tests for teloclip.analysis module.

Tests for overhang statistics collection, outlier detection,
and homopolymer analysis functionality.
"""

import pytest

from teloclip.analysis import (
    ContigStats,
    OverhangInfo,
    calculate_overhang_statistics,
    collect_overhang_stats,
    detect_homopolymer_runs,
    identify_outlier_contigs,
    rank_overhangs_by_length,
    select_best_overhang,
)


class TestOverhangInfo:
    """Test the OverhangInfo dataclass."""

    def test_overhang_info_creation(self):
        """Test basic OverhangInfo creation."""
        overhang = OverhangInfo(
            sequence='ATCG',
            length=4,
            alignment_pos=1,
            alignment_end=100,
            read_name='read1',
            is_left=True,
            clip_length=4,
            anchor_length=96,
        )

        assert overhang.sequence == 'ATCG'
        assert overhang.length == 4
        assert overhang.is_left is True
        assert overhang.read_name == 'read1'


class TestContigStats:
    """Test the ContigStats dataclass."""

    def test_contig_stats_creation(self):
        """Test basic ContigStats creation."""
        stats = ContigStats('contig1', 1000)

        assert stats.contig_name == 'contig1'
        assert stats.contig_length == 1000
        assert stats.left_count == 0
        assert stats.right_count == 0
        assert stats.left_total_length == 0
        assert stats.right_total_length == 0

    def test_contig_stats_properties(self):
        """Test ContigStats computed properties."""
        overhang1 = OverhangInfo('ATCG', 4, 1, 100, 'read1', True, 4, 96)
        overhang2 = OverhangInfo('GCTA', 4, 1, 100, 'read2', True, 4, 96)
        overhang3 = OverhangInfo('TTAA', 4, 995, 1000, 'read3', False, 4, 96)

        stats = ContigStats('contig1', 1000)
        stats.left_overhangs = [overhang1, overhang2]
        stats.right_overhangs = [overhang3]

        assert stats.left_count == 2
        assert stats.right_count == 1
        assert stats.left_total_length == 8
        assert stats.right_total_length == 4


class TestCollectOverhangStats:
    """Test the collect_overhang_stats function."""

    def test_collect_basic_stats(self):
        """Test basic overhang statistics collection."""
        sam_lines = [
            '@HD\tVN:1.0\tSO:coordinate',
            'read1\t0\tcontig1\t1\t60\t20S80M\t*\t0\t0\tAAAAAAAAAAAAAAAAAATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t*',
            'read2\t0\tcontig1\t990\t60\t80M20S\t*\t0\t0\tTCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGGGGGGGGGGGGGGGGGGGGG\t*',
        ]

        contig_dict = {'contig1': 1000}

        stats = collect_overhang_stats(iter(sam_lines), contig_dict)

        assert 'contig1' in stats
        assert stats['contig1'].left_count == 1
        assert stats['contig1'].right_count == 1
        assert stats['contig1'].left_overhangs[0].sequence == 'AAAAAAAAAAAAAAAAAATC'
        assert stats['contig1'].right_overhangs[0].sequence == 'GGGGGGGGGGGGGGGGGGGG'

    def test_collect_stats_skip_unmapped(self):
        """Test that unmapped reads are skipped."""
        sam_lines = [
            'read1\t4\t*\t0\t0\t*\t*\t0\t0\tATCGATCGATCGATCG\t*',  # unmapped
            'read2\t0\tcontig1\t1\t60\t20S80M\t*\t0\t0\tAAAAAAAAAAAAAAAAAATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t*',
        ]

        contig_dict = {'contig1': 1000}

        stats = collect_overhang_stats(iter(sam_lines), contig_dict)

        assert stats['contig1'].left_count == 1

    def test_collect_stats_skip_secondary(self):
        """Test that secondary alignments are skipped."""
        sam_lines = [
            'read1\t256\tcontig1\t1\t60\t20S80M\t*\t0\t0\tAAAAAAAAAAAAAAAAAATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t*',  # secondary
            'read2\t0\tcontig1\t1\t60\t20S80M\t*\t0\t0\tAAAAAAAAAAAAAAAAAATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t*',
        ]

        contig_dict = {'contig1': 1000}

        stats = collect_overhang_stats(iter(sam_lines), contig_dict)

        assert stats['contig1'].left_count == 1


class TestCalculateOverhangStatistics:
    """Test the calculate_overhang_statistics function."""

    def test_calculate_statistics_empty(self):
        """Test statistics calculation with empty data."""
        stats_dict = {'contig1': ContigStats('contig1', 1000)}

        result = calculate_overhang_statistics(stats_dict)

        assert result['left']['mean'] == 0.0
        assert result['right']['mean'] == 0.0
        assert result['combined']['mean'] == 0.0

    def test_calculate_statistics_with_data(self):
        """Test statistics calculation with actual data."""
        overhang1 = OverhangInfo('ATCG', 4, 1, 100, 'read1', True, 4, 96)
        overhang2 = OverhangInfo('ATCGAA', 6, 1, 100, 'read2', True, 6, 94)
        overhang3 = OverhangInfo('TTAA', 4, 995, 1000, 'read3', False, 4, 96)

        stats = ContigStats('contig1', 1000)
        stats.left_overhangs = [overhang1, overhang2]
        stats.right_overhangs = [overhang3]

        stats_dict = {'contig1': stats}

        result = calculate_overhang_statistics(stats_dict)

        assert result['left']['mean'] == 5.0  # (4 + 6) / 2
        assert result['right']['mean'] == 4.0
        assert result['combined']['mean'] == pytest.approx(
            4.67, rel=1e-2
        )  # (4 + 6 + 4) / 3


class TestIdentifyOutlierContigs:
    """Test the identify_outlier_contigs function."""

    def test_identify_outliers_no_outliers(self):
        """Test outlier detection with no outliers."""
        # Create similar overhang counts
        stats_dict = {}
        for i in range(5):
            stats = ContigStats(f'contig{i}', 1000)
            # Add 3 overhangs to each side for each contig
            for j in range(3):
                left_oh = OverhangInfo('ATCG', 4, 1, 100, f'read{j}', True, 4, 96)
                right_oh = OverhangInfo(
                    'TTAA', 4, 995, 1000, f'read{j + 3}', False, 4, 96
                )
                stats.left_overhangs.append(left_oh)
                stats.right_overhangs.append(right_oh)
            stats_dict[f'contig{i}'] = stats

        outliers = identify_outlier_contigs(stats_dict, threshold=2.0)

        assert len(outliers['left_outliers']) == 0
        assert len(outliers['right_outliers']) == 0

    def test_identify_outliers_with_outliers(self):
        """Test outlier detection with actual outliers."""
        stats_dict = {}

        # Create normal contigs with 3 overhangs each
        for i in range(4):
            stats = ContigStats(f'contig{i}', 1000)
            for j in range(3):
                left_oh = OverhangInfo('ATCG', 4, 1, 100, f'read{j}', True, 4, 96)
                right_oh = OverhangInfo(
                    'TTAA', 4, 995, 1000, f'read{j + 3}', False, 4, 96
                )
                stats.left_overhangs.append(left_oh)
                stats.right_overhangs.append(right_oh)
            stats_dict[f'contig{i}'] = stats

        # Create outlier contig with many overhangs
        outlier_stats = ContigStats('outlier', 1000)
        for j in range(30):  # Much higher than others - increased from 20 to 30
            left_oh = OverhangInfo('ATCG', 4, 1, 100, f'read{j}', True, 4, 96)
            outlier_stats.left_overhangs.append(left_oh)
        stats_dict['outlier'] = outlier_stats

        outliers = identify_outlier_contigs(stats_dict, threshold=1.5)

        assert 'outlier' in outliers['left_outliers']


class TestRankOverhangsByLength:
    """Test the rank_overhangs_by_length function."""

    def test_rank_overhangs(self):
        """Test ranking overhangs by length."""
        overhangs = [
            OverhangInfo('AT', 2, 1, 100, 'read1', True, 2, 98),
            OverhangInfo('ATCGATCG', 8, 1, 100, 'read2', True, 8, 92),
            OverhangInfo('ATCG', 4, 1, 100, 'read3', True, 4, 96),
        ]

        ranked = rank_overhangs_by_length(overhangs)

        assert ranked[0].length == 8
        assert ranked[1].length == 4
        assert ranked[2].length == 2

    def test_rank_empty_list(self):
        """Test ranking empty list."""
        ranked = rank_overhangs_by_length([])
        assert len(ranked) == 0


class TestDetectHomopolymerRuns:
    """Test the detect_homopolymer_runs function."""

    def test_detect_homopolymer_basic(self):
        """Test basic homopolymer detection."""
        sequence = 'ATCG' + 'A' * 60 + 'TCGA'

        runs = detect_homopolymer_runs(sequence, min_length=50)

        assert len(runs) == 1
        assert runs[0][0] == 'A'  # nucleotide
        assert runs[0][3] == 60  # length

    def test_detect_multiple_homopolymers(self):
        """Test detection of multiple homopolymer runs."""
        sequence = 'A' * 60 + 'TCGA' + 'G' * 55 + 'ATCG'

        runs = detect_homopolymer_runs(sequence, min_length=50)

        assert len(runs) == 2
        assert runs[0][0] == 'A'
        assert runs[0][3] == 60
        assert runs[1][0] == 'G'
        assert runs[1][3] == 55

    def test_detect_no_homopolymers(self):
        """Test sequence with no long homopolymer runs."""
        sequence = 'ATCGATCGATCGATCG' * 10

        runs = detect_homopolymer_runs(sequence, min_length=50)

        assert len(runs) == 0

    def test_detect_empty_sequence(self):
        """Test empty sequence."""
        runs = detect_homopolymer_runs('', min_length=50)
        assert len(runs) == 0


class TestSelectBestOverhang:
    """Test the select_best_overhang function."""

    def test_select_longest_overhang(self):
        """Test selection of longest overhang."""
        overhangs = [
            OverhangInfo('AT', 2, 1, 100, 'read1', True, 2, 98),
            OverhangInfo('ATCGATCG', 8, 1, 100, 'read2', True, 8, 92),
            OverhangInfo('ATCG', 4, 1, 100, 'read3', True, 4, 96),
        ]

        best = select_best_overhang(overhangs, min_extension=1, max_homopolymer=50)

        assert best.length == 8
        assert best.read_name == 'read2'

    def test_select_with_min_extension_filter(self):
        """Test selection with minimum extension requirement."""
        overhangs = [
            OverhangInfo('AT', 2, 1, 100, 'read1', True, 2, 98),
            OverhangInfo('ATCGATCG', 8, 1, 100, 'read2', True, 8, 92),
            OverhangInfo('ATCG', 4, 1, 100, 'read3', True, 4, 96),
        ]

        best = select_best_overhang(overhangs, min_extension=5, max_homopolymer=50)

        assert best.length == 8
        assert best.read_name == 'read2'

    def test_select_with_homopolymer_filtering(self):
        """Test selection avoiding homopolymer runs."""
        overhangs = [
            OverhangInfo(
                'A' * 60, 60, 1, 100, 'read1', True, 60, 40
            ),  # Has homopolymer
            OverhangInfo('ATCGATCG', 8, 1, 100, 'read2', True, 8, 92),  # Clean
            OverhangInfo('ATCG', 4, 1, 100, 'read3', True, 4, 96),  # Clean but shorter
        ]

        best = select_best_overhang(overhangs, min_extension=1, max_homopolymer=50)

        # Should select the clean sequence even though it's shorter
        assert best.length == 8
        assert best.read_name == 'read2'

    def test_select_empty_list(self):
        """Test selection from empty list."""
        best = select_best_overhang([], min_extension=1, max_homopolymer=50)
        assert best is None

    def test_select_no_candidates_meet_criteria(self):
        """Test when no overhangs meet the criteria."""
        overhangs = [
            OverhangInfo('AT', 2, 1, 100, 'read1', True, 2, 98),
            OverhangInfo('A', 1, 1, 100, 'read2', True, 1, 99),
        ]

        best = select_best_overhang(overhangs, min_extension=5, max_homopolymer=50)
        assert best is None
