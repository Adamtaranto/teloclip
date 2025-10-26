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
            contig_name='test_contig',
        )

        assert overhang.sequence == 'ATCG'
        assert overhang.length == 4
        assert overhang.is_left is True
        assert overhang.read_name == 'read1'
        assert overhang.contig_name == 'test_contig'


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
        overhang1 = OverhangInfo('ATCG', 4, 1, 100, 'read1', True, 4, 96, 'contig1')
        overhang2 = OverhangInfo('GCTA', 4, 1, 100, 'read2', True, 4, 96, 'contig1')
        overhang3 = OverhangInfo('TTAA', 4, 995, 1000, 'read3', False, 4, 96, 'contig1')

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
        overhang1 = OverhangInfo('ATCG', 4, 1, 100, 'read1', True, 4, 96, 'contig1')
        overhang2 = OverhangInfo('ATCGAA', 6, 1, 100, 'read2', True, 6, 94, 'contig1')
        overhang3 = OverhangInfo('TTAA', 4, 995, 1000, 'read3', False, 4, 96, 'contig1')

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
                left_oh = OverhangInfo(
                    'ATCG', 4, 1, 100, f'read{j}', True, 4, 96, f'contig{i}'
                )
                right_oh = OverhangInfo(
                    'TTAA', 4, 995, 1000, f'read{j + 3}', False, 4, 96, f'contig{i}'
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
                left_oh = OverhangInfo(
                    'ATCG', 4, 1, 100, f'read{j}', True, 4, 96, f'contig{i}'
                )
                right_oh = OverhangInfo(
                    'TTAA', 4, 995, 1000, f'read{j + 3}', False, 4, 96, f'contig{i}'
                )
                stats.left_overhangs.append(left_oh)
                stats.right_overhangs.append(right_oh)
            stats_dict[f'contig{i}'] = stats

        # Create outlier contig with many overhangs
        outlier_stats = ContigStats('outlier', 1000)
        for j in range(30):  # Much higher than others - increased from 20 to 30
            left_oh = OverhangInfo(
                'ATCG', 4, 1, 100, f'read{j}', True, 4, 96, 'outlier'
            )
            outlier_stats.left_overhangs.append(left_oh)
        stats_dict['outlier'] = outlier_stats

        outliers = identify_outlier_contigs(stats_dict, threshold=1.5)

        assert 'outlier' in outliers['left_outliers']


class TestRankOverhangsByLength:
    """Test the rank_overhangs_by_length function."""

    def test_rank_overhangs(self):
        """Test ranking overhangs by length."""
        overhangs = [
            OverhangInfo('AT', 2, 1, 100, 'read1', True, 2, 98, 'test_contig'),
            OverhangInfo('ATCGATCG', 8, 1, 100, 'read2', True, 8, 92, 'test_contig'),
            OverhangInfo('ATCG', 4, 1, 100, 'read3', True, 4, 96, 'test_contig'),
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
            OverhangInfo('AT', 2, 1, 100, 'read1', True, 2, 98, 'test_contig'),
            OverhangInfo('ATCGATCG', 8, 1, 100, 'read2', True, 8, 92, 'test_contig'),
            OverhangInfo('ATCG', 4, 1, 100, 'read3', True, 4, 96, 'test_contig'),
        ]

        best = select_best_overhang(overhangs, min_extension=1, max_homopolymer=50)

        assert best.length == 8
        assert best.read_name == 'read2'

    def test_select_with_min_extension_filter(self):
        """Test selection with minimum extension requirement."""
        overhangs = [
            OverhangInfo('AT', 2, 1, 100, 'read1', True, 2, 98, 'test_contig'),
            OverhangInfo('ATCGATCG', 8, 1, 100, 'read2', True, 8, 92, 'test_contig'),
            OverhangInfo('ATCG', 4, 1, 100, 'read3', True, 4, 96, 'test_contig'),
        ]

        best = select_best_overhang(overhangs, min_extension=5, max_homopolymer=50)

        assert best.length == 8
        assert best.read_name == 'read2'

    def test_select_with_homopolymer_filtering(self):
        """Test selection avoiding homopolymer runs."""
        overhangs = [
            OverhangInfo(
                'A' * 60, 60, 1, 100, 'read1', True, 60, 40, 'test_contig'
            ),  # Has homopolymer
            OverhangInfo(
                'ATCGATCG', 8, 1, 100, 'read2', True, 8, 92, 'test_contig'
            ),  # Clean
            OverhangInfo(
                'ATCG', 4, 1, 100, 'read3', True, 4, 96, 'test_contig'
            ),  # Clean but shorter
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
            OverhangInfo('AT', 2, 1, 100, 'read1', True, 2, 98, 'test_contig'),
            OverhangInfo('A', 1, 1, 100, 'read2', True, 1, 99, 'test_contig'),
        ]

        best = select_best_overhang(overhangs, min_extension=5, max_homopolymer=50)
        assert best is None


class TestDualEndOverhang:
    """Test cases for reads that overhang both ends of a contig."""

    def setup_method(self):
        """Set up common test data."""
        from teloclip.streaming_analysis import process_single_contig_extension

        self.process_extension = process_single_contig_extension

    def create_dual_end_contig_stats(self, read_name: str = 'dual_overhang_read'):
        """Create ContigStats with a single read overhanging both ends."""
        contig_stats = ContigStats('short_contig', 500)  # Short contig

        # Single read that overhangs both left and right ends
        # Left overhang: read starts before contig, extends into contig
        left_overhang = OverhangInfo(
            sequence='TTAGGGTTAGGGTTAGGG',  # 18bp left overhang
            length=18,
            alignment_pos=1,  # Alignment starts at contig position 1
            alignment_end=450,  # Alignment extends to position 450
            read_name=read_name,
            is_left=True,
            clip_length=18,
            anchor_length=450,  # 450bp anchored in contig
            contig_name='short_contig',
        )

        # Right overhang: same read continues past contig end
        right_overhang = OverhangInfo(
            sequence='CCCTAACCCTAACCCTAA',  # 18bp right overhang
            length=18,
            alignment_pos=51,  # Alignment starts at position 51
            alignment_end=500,  # Alignment extends to contig end
            read_name=read_name,  # Same read name!
            is_left=False,
            clip_length=18,
            anchor_length=450,  # 450bp anchored in contig
            contig_name='short_contig',
        )

        contig_stats.left_overhangs = [left_overhang]
        contig_stats.right_overhangs = [right_overhang]

        return contig_stats

    def test_dual_end_overhang_same_read(self):
        """Test extension when the same read overhangs both ends."""
        contig_stats = self.create_dual_end_contig_stats('spanning_read_001')
        original_sequence = 'A' * 500  # Simple 500bp sequence

        result = self.process_extension(
            contig_name='short_contig',
            contig_stats=contig_stats,
            original_sequence=original_sequence,
            min_extension=1,
            max_homopolymer=100,
            motif_patterns=None,
            dry_run=False,
        )

        assert result is not None, 'Should successfully extend with dual-end overhang'

        # Verify both extensions are recorded
        ext_info = result.extension_info
        assert ext_info.get('has_left_extension', False), 'Should have left extension'
        assert ext_info.get('has_right_extension', False), 'Should have right extension'

        # Verify extension details
        assert ext_info.get('left_overhang_length') == 18, (
            'Left extension should be 18bp'
        )
        assert ext_info.get('right_overhang_length') == 18, (
            'Right extension should be 18bp'
        )
        assert ext_info.get('left_read_name') == 'spanning_read_001', (
            'Left read name should match'
        )
        assert ext_info.get('right_read_name') == 'spanning_read_001', (
            'Right read name should match'
        )

        # Verify final sequence length
        expected_length = 500 + 18 + 18  # Original + left + right
        assert result.extension_info['final_length'] == expected_length
        assert len(result.extended_sequence) == expected_length

        # Verify sequence structure
        assert result.extended_sequence.startswith('TTAGGGTTAGGGTTAGGG'), (
            'Should start with left overhang'
        )
        assert result.extended_sequence.endswith('CCCTAACCCTAACCCTAA'), (
            'Should end with right overhang'
        )

    def test_dual_end_with_competing_overhangs(self):
        """Test when the dual-end read competes with other shorter overhangs."""
        contig_stats = self.create_dual_end_contig_stats('long_spanning_read')

        # Add competing shorter overhangs
        short_left = OverhangInfo(
            sequence='TTAG',  # 4bp - shorter
            length=4,
            alignment_pos=1,
            alignment_end=200,
            read_name='short_left_read',
            is_left=True,
            clip_length=4,
            anchor_length=199,
            contig_name='test_contig',
        )

        short_right = OverhangInfo(
            sequence='CCTA',  # 4bp - shorter
            length=4,
            alignment_pos=301,
            alignment_end=500,
            read_name='short_right_read',
            is_left=False,
            clip_length=4,
            anchor_length=199,
            contig_name='test_contig',
        )

        contig_stats.left_overhangs.append(short_left)
        contig_stats.right_overhangs.append(short_right)

        original_sequence = 'A' * 500

        result = self.process_extension(
            contig_name='short_contig',
            contig_stats=contig_stats,
            original_sequence=original_sequence,
            min_extension=1,
            max_homopolymer=100,
            motif_patterns=None,
            dry_run=False,
        )

        assert result is not None

        # The longer overhangs from the spanning read should be selected
        ext_info = result.extension_info
        assert ext_info.get('left_overhang_length') == 18, (
            'Should select longer left overhang'
        )
        assert ext_info.get('right_overhang_length') == 18, (
            'Should select longer right overhang'
        )
        assert ext_info.get('left_read_name') == 'long_spanning_read', (
            'Should use spanning read for left'
        )
        assert ext_info.get('right_read_name') == 'long_spanning_read', (
            'Should use spanning read for right'
        )

    def test_dual_end_dry_run(self):
        """Test dual-end extension in dry run mode."""
        contig_stats = self.create_dual_end_contig_stats('dry_run_read')
        original_sequence = 'A' * 500

        result = self.process_extension(
            contig_name='short_contig',
            contig_stats=contig_stats,
            original_sequence=original_sequence,
            min_extension=1,
            max_homopolymer=100,
            motif_patterns=None,
            dry_run=True,
        )

        assert result is not None

        # In dry run, sequence should be unchanged but extensions should be recorded
        assert result.extended_sequence == original_sequence, (
            'Sequence should be unchanged in dry run'
        )

        ext_info = result.extension_info
        assert ext_info.get('has_left_extension', False), 'Should plan left extension'
        assert ext_info.get('has_right_extension', False), 'Should plan right extension'
        assert ext_info.get('final_length') == 536, (
            'Should calculate final length correctly'
        )

    def test_dual_end_with_homopolymer_filtering(self):
        """Test dual-end extension when overhangs contain homopolymer runs."""
        contig_stats = ContigStats('short_contig', 500)

        # Create overhangs with homopolymer runs that exceed threshold
        left_with_homopolymer = OverhangInfo(
            sequence='A' * 60,  # Long homopolymer run
            length=60,
            alignment_pos=1,
            alignment_end=400,
            read_name='homopolymer_read',
            is_left=True,
            clip_length=60,
            anchor_length=399,
            contig_name='test_contig',
        )

        right_with_homopolymer = OverhangInfo(
            sequence='T' * 60,  # Long homopolymer run
            length=60,
            alignment_pos=101,
            alignment_end=500,
            read_name='homopolymer_read',
            is_left=False,
            clip_length=60,
            anchor_length=399,
            contig_name='test_contig',
        )

        contig_stats.left_overhangs = [left_with_homopolymer]
        contig_stats.right_overhangs = [right_with_homopolymer]

        original_sequence = 'A' * 500

        result = self.process_extension(
            contig_name='short_contig',
            contig_stats=contig_stats,
            original_sequence=original_sequence,
            min_extension=1,
            max_homopolymer=50,  # Threshold lower than overhang homopolymer length
            motif_patterns=None,
            dry_run=False,
        )

        # Should not extend due to homopolymer runs exceeding threshold
        assert result is None, (
            'Should not extend when all overhangs contain homopolymers exceeding threshold'
        )

    def test_dual_end_minimum_extension_filtering(self):
        """Test dual-end extension with minimum extension length filtering."""
        contig_stats = ContigStats('short_contig', 500)

        # Create short overhangs that don't meet minimum extension requirement
        short_left = OverhangInfo(
            sequence='AT',  # 2bp - below minimum
            length=2,
            alignment_pos=1,
            alignment_end=400,
            read_name='short_read',
            is_left=True,
            clip_length=2,
            anchor_length=399,
            contig_name='test_contig',
        )

        short_right = OverhangInfo(
            sequence='GC',  # 2bp - below minimum
            length=2,
            alignment_pos=101,
            alignment_end=500,
            read_name='short_read',
            is_left=False,
            clip_length=2,
            anchor_length=399,
            contig_name='test_contig',
        )

        contig_stats.left_overhangs = [short_left]
        contig_stats.right_overhangs = [short_right]

        original_sequence = 'A' * 500

        result = self.process_extension(
            contig_name='short_contig',
            contig_stats=contig_stats,
            original_sequence=original_sequence,
            min_extension=5,  # Minimum higher than overhang lengths
            max_homopolymer=100,
            motif_patterns=None,
            dry_run=False,
        )

        # Should return None because no overhangs meet minimum extension requirement
        assert result is None, (
            'Should not extend when overhangs are below minimum length'
        )

    def test_dual_end_asymmetric_filtering(self):
        """Test when only one end meets extension criteria."""
        contig_stats = ContigStats('short_contig', 500)

        # Good left overhang
        good_left = OverhangInfo(
            sequence='TTAGGGTTAGGG',  # 12bp - meets criteria
            length=12,
            alignment_pos=1,
            alignment_end=400,
            read_name='good_read',
            is_left=True,
            clip_length=12,
            anchor_length=399,
            contig_name='test_contig',
        )

        # Bad right overhang (homopolymer)
        bad_right = OverhangInfo(
            sequence='A' * 60,  # Homopolymer run
            length=60,
            alignment_pos=101,
            alignment_end=500,
            read_name='bad_read',
            is_left=False,
            clip_length=60,
            anchor_length=399,
            contig_name='test_contig',
        )

        contig_stats.left_overhangs = [good_left]
        contig_stats.right_overhangs = [bad_right]

        original_sequence = 'A' * 500

        result = self.process_extension(
            contig_name='short_contig',
            contig_stats=contig_stats,
            original_sequence=original_sequence,
            min_extension=1,
            max_homopolymer=50,  # Will filter out right overhang
            motif_patterns=None,
            dry_run=False,
        )

        assert result is not None, 'Should still extend left end only'

        ext_info = result.extension_info
        assert ext_info.get('has_left_extension', False), 'Should have left extension'
        assert not ext_info.get('has_right_extension', False), (
            'Should NOT have right extension due to homopolymer filtering'
        )

        # Should only extend left end (clean overhang)
        expected_length = 500 + 12  # Original + left only (right rejected)
        assert result.extension_info['final_length'] == expected_length
        assert len(result.extended_sequence) == expected_length

    def test_dual_end_single_read_different_qualities(self):
        """Test dual-end overhang from single read with different anchor qualities."""
        contig_stats = ContigStats('short_contig', 300)  # Even shorter contig

        # Same read with different anchor qualities on each end
        # This simulates a read that spans most of a short contig
        left_overhang = OverhangInfo(
            sequence='TTAGGGTTAGGGTTAGGGTTAGGG',  # 24bp left overhang
            length=24,
            alignment_pos=1,  # Perfect start alignment
            alignment_end=250,  # Good coverage
            read_name='spanning_read_high_quality',
            is_left=True,
            clip_length=24,
            anchor_length=249,  # High anchor quality
            contig_name='short_contig',
        )

        right_overhang = OverhangInfo(
            sequence='CCCTAACCCTAACCCTAACCCTAA',  # 24bp right overhang
            length=24,
            alignment_pos=51,  # Starts later in contig
            alignment_end=300,  # Perfect end alignment
            read_name='spanning_read_high_quality',  # Same read
            is_left=False,
            clip_length=24,
            anchor_length=249,  # High anchor quality
            contig_name='short_contig',
        )

        contig_stats.left_overhangs = [left_overhang]
        contig_stats.right_overhangs = [right_overhang]

        original_sequence = 'G' * 300  # 300bp sequence different from test above

        result = self.process_extension(
            contig_name='short_contig',
            contig_stats=contig_stats,
            original_sequence=original_sequence,
            min_extension=1,
            max_homopolymer=100,
            motif_patterns=None,
            dry_run=False,
        )

        assert result is not None, (
            'Should successfully extend with high-quality dual overhang'
        )

        # Both extensions should be recorded
        ext_info = result.extension_info
        assert ext_info.get('has_left_extension', False), 'Should have left extension'
        assert ext_info.get('has_right_extension', False), 'Should have right extension'

        # Same read should be used for both extensions
        assert ext_info.get('left_read_name') == 'spanning_read_high_quality'
        assert ext_info.get('right_read_name') == 'spanning_read_high_quality'

        # Both extensions should be 24bp
        assert ext_info.get('left_overhang_length') == 24
        assert ext_info.get('right_overhang_length') == 24

        # Final sequence should be original + both overhangs
        expected_length = 300 + 24 + 24  # 348bp total
        assert result.extension_info['final_length'] == expected_length
        assert len(result.extended_sequence) == expected_length

        # Verify sequence structure - should have overhangs on both ends
        assert result.extended_sequence.startswith('TTAGGGTTAGGGTTAGGGTTAGGG'), (
            'Should start with left overhang'
        )
        assert result.extended_sequence.endswith('CCCTAACCCTAACCCTAACCCTAA'), (
            'Should end with right overhang'
        )

        # Middle should be the original sequence
        middle_seq = result.extended_sequence[24:-24]  # Remove overhangs
        assert middle_seq == original_sequence, 'Middle should be original sequence'

    def test_single_distinct_left_overhang_only(self):
        """Test when only left end has a distinct best overhang (no right overhangs)."""
        contig_stats = ContigStats('left_only_contig', 800)

        # Single best left overhang
        best_left = OverhangInfo(
            sequence='TTAGGGTTAGGGTTAGGGTTAGGG',  # 24bp distinctive sequence
            length=24,
            alignment_pos=1,
            alignment_end=600,
            read_name='distinct_left_read',
            is_left=True,
            clip_length=24,
            anchor_length=599,
            contig_name='left_only_contig',
        )

        contig_stats.left_overhangs = [best_left]
        contig_stats.right_overhangs = []  # No right overhangs

        original_sequence = 'C' * 800

        result = self.process_extension(
            contig_name='left_only_contig',
            contig_stats=contig_stats,
            original_sequence=original_sequence,
            min_extension=1,
            max_homopolymer=100,
            motif_patterns=None,
            dry_run=False,
        )

        assert result is not None, 'Should extend left end even with no right overhangs'

        ext_info = result.extension_info
        assert ext_info.get('has_left_extension', False), 'Should have left extension'
        assert not ext_info.get('has_right_extension', False), (
            'Should NOT have right extension'
        )

        # Verify left extension details
        assert ext_info.get('left_overhang_length') == 24
        assert ext_info.get('left_read_name') == 'distinct_left_read'

        # Final sequence should be original + left extension only
        expected_length = 800 + 24  # 824bp total
        assert result.extension_info['final_length'] == expected_length
        assert len(result.extended_sequence) == expected_length

        # Should start with left overhang, end with original sequence
        assert result.extended_sequence.startswith('TTAGGGTTAGGGTTAGGGTTAGGG')
        assert result.extended_sequence.endswith(
            'C' * 100
        )  # Last 100bp should be original

    def test_single_distinct_right_overhang_only(self):
        """Test when only right end has a distinct best overhang (no left overhangs)."""
        contig_stats = ContigStats('right_only_contig', 800)

        # Single best right overhang
        best_right = OverhangInfo(
            sequence='CCCTAACCCTAACCCTAACCCTAA',  # 24bp distinctive sequence
            length=24,
            alignment_pos=201,
            alignment_end=800,
            read_name='distinct_right_read',
            is_left=False,
            clip_length=24,
            anchor_length=599,
            contig_name='right_only_contig',
        )

        contig_stats.left_overhangs = []  # No left overhangs
        contig_stats.right_overhangs = [best_right]

        original_sequence = 'G' * 800

        result = self.process_extension(
            contig_name='right_only_contig',
            contig_stats=contig_stats,
            original_sequence=original_sequence,
            min_extension=1,
            max_homopolymer=100,
            motif_patterns=None,
            dry_run=False,
        )

        assert result is not None, 'Should extend right end even with no left overhangs'

        ext_info = result.extension_info
        assert not ext_info.get('has_left_extension', False), (
            'Should NOT have left extension'
        )
        assert ext_info.get('has_right_extension', False), 'Should have right extension'

        # Verify right extension details
        assert ext_info.get('right_overhang_length') == 24
        assert ext_info.get('right_read_name') == 'distinct_right_read'

        # Final sequence should be original + right extension only
        expected_length = 800 + 24  # 824bp total
        assert result.extension_info['final_length'] == expected_length
        assert len(result.extended_sequence) == expected_length

        # Should start with original sequence, end with right overhang
        assert result.extended_sequence.startswith(
            'G' * 100
        )  # First 100bp should be original
        assert result.extended_sequence.endswith('CCCTAACCCTAACCCTAACCCTAA')

    def test_single_read_longest_on_both_ends(self):
        """Test single read spanning contig that's definitively the longest on BOTH ends."""
        contig_stats = ContigStats('spanning_contig', 400)

        # Single spanning read with LONG overhangs on both ends
        spanning_left = OverhangInfo(
            sequence='TTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGG',  # 36bp - LONGEST
            length=36,
            alignment_pos=1,
            alignment_end=350,
            read_name='longest_spanning_read',
            is_left=True,
            clip_length=36,
            anchor_length=349,
            contig_name='spanning_contig',
        )

        spanning_right = OverhangInfo(
            sequence='CCCTAACCCTAACCCTAACCCTAACCCTAACCCTAA',  # 36bp - LONGEST
            length=36,
            alignment_pos=51,
            alignment_end=400,
            read_name='longest_spanning_read',  # Same read!
            is_left=False,
            clip_length=36,
            anchor_length=349,
            contig_name='spanning_contig',
        )

        # Add multiple competing shorter overhangs to test selection
        competitors = []
        for i in range(5):
            # Shorter left competitors
            left_competitor = OverhangInfo(
                sequence='TTAG' * (2 + i),  # 8-20bp sequences
                length=8 + 4 * i,
                alignment_pos=1,
                alignment_end=200 + i * 10,
                read_name=f'short_left_{i}',
                is_left=True,
                clip_length=8 + 4 * i,
                anchor_length=199 + i * 10,
                contig_name='spanning_contig',
            )

            # Shorter right competitors
            right_competitor = OverhangInfo(
                sequence='CCTA' * (2 + i),  # 8-20bp sequences
                length=8 + 4 * i,
                alignment_pos=201 + i * 10,
                alignment_end=400,
                read_name=f'short_right_{i}',
                is_left=False,
                clip_length=8 + 4 * i,
                anchor_length=199 + i * 10,
                contig_name='spanning_contig',
            )

            competitors.extend([left_competitor, right_competitor])

        # Set up contig stats with spanning read + competitors
        contig_stats.left_overhangs = [spanning_left] + [
            c for c in competitors if c.is_left
        ]
        contig_stats.right_overhangs = [spanning_right] + [
            c for c in competitors if not c.is_left
        ]

        original_sequence = 'A' * 400

        result = self.process_extension(
            contig_name='spanning_contig',
            contig_stats=contig_stats,
            original_sequence=original_sequence,
            min_extension=1,
            max_homopolymer=100,
            motif_patterns=None,
            dry_run=False,
        )

        assert result is not None, (
            'Should successfully extend with longest spanning read'
        )

        # Verify both extensions use the same spanning read (longest on both ends)
        ext_info = result.extension_info
        assert ext_info.get('has_left_extension', False), 'Should have left extension'
        assert ext_info.get('has_right_extension', False), 'Should have right extension'

        # Verify the longest overhangs were selected
        assert ext_info.get('left_overhang_length') == 36, (
            'Should select longest left overhang'
        )
        assert ext_info.get('right_overhang_length') == 36, (
            'Should select longest right overhang'
        )

        # Verify same read used for both (the spanning read)
        assert ext_info.get('left_read_name') == 'longest_spanning_read'
        assert ext_info.get('right_read_name') == 'longest_spanning_read'

        # Verify final sequence incorporates both long extensions
        expected_length = 400 + 36 + 36  # 472bp total
        assert result.extension_info['final_length'] == expected_length
        assert len(result.extended_sequence) == expected_length

        # Verify sequence structure with longest overhangs
        assert result.extended_sequence.startswith(
            'TTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGG'
        ), 'Should start with longest left overhang'
        assert result.extended_sequence.endswith(
            'CCCTAACCCTAACCCTAACCCTAACCCTAACCCTAA'
        ), 'Should end with longest right overhang'
