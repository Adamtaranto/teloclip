"""
Unit tests for teloclip.core.analysis module.

Tests for overhang statistics collection, outlier detection,
and homopolymer analysis functionality.
"""

import pytest

from teloclip.core.analysis import (
    ContigStats,
    OverhangInfo,
    calculate_overhang_statistics,
    detect_homopolymer_runs,
    flag_anomalous_overhang_coverage,
    rank_overhangs_by_gain,
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


def make_stats(counts):
    """
    Build a stats mapping from per-contig (left, right) overhang counts.

    Parameters
    ----------
    counts : dict
        Mapping of contig name to a ``(left_count, right_count)`` tuple.

    Returns
    -------
    dict
        Mapping of contig name to a populated ContigStats.
    """
    stats_dict = {}
    for name, (left, right) in counts.items():
        stats = ContigStats(name, 1000)
        stats.left_overhangs = [
            OverhangInfo('ATCG', 4, 1, 100, f'l{i}', True, 4, 96, name, 4)
            for i in range(left)
        ]
        stats.right_overhangs = [
            OverhangInfo('TTAA', 4, 995, 1000, f'r{i}', False, 4, 96, name, 4)
            for i in range(right)
        ]
        stats_dict[name] = stats
    return stats_dict


class TestFlagAnomalousOverhangCoverage:
    """Test flagging of contigs with unusual overhang coverage."""

    def test_uniform_coverage_flags_nothing(self):
        """Test that an assembly with even coverage produces no flags."""
        stats = make_stats({f'contig{i}': (3, 3) for i in range(10)})

        result = flag_anomalous_overhang_coverage(stats)

        assert result['left_outliers'] == []
        assert result['right_outliers'] == []

    def test_high_tail_is_flagged(self):
        """Test that a contig hoarding clipped reads is flagged."""
        counts = {f'contig{i}': (3, 3) for i in range(10)}
        counts['rdna_array'] = (400, 3)
        stats = make_stats(counts)

        result = flag_anomalous_overhang_coverage(stats)

        assert result['left_outliers'] == ['rdna_array']
        assert result['right_outliers'] == []

    def test_low_tail_is_not_flagged(self):
        """Test that a contig with unusually few overhangs is not flagged.

        Sparse evidence is not an anomaly. The previous two-tailed z-score
        excluded these contigs from extension identically to the high tail,
        which is the opposite of useful.
        """
        counts = {f'contig{i}': (40, 40) for i in range(10)}
        counts['quiet_contig'] = (1, 1)
        stats = make_stats(counts)

        result = flag_anomalous_overhang_coverage(stats)

        assert 'quiet_contig' not in result['left_outliers']
        assert 'quiet_contig' not in result['right_outliers']

    def test_declines_to_judge_a_small_assembly(self):
        """Test that too few contigs produces no flags at all.

        With a handful of contigs the spread carries too little information.
        The previous mean/stdev implementation could not exceed a z-score of 2
        with four contigs regardless of how extreme one of them was, so the
        feature silently did nothing on small assemblies.
        """
        counts = {f'contig{i}': (3, 3) for i in range(4)}
        counts['extreme'] = (10_000, 3)
        stats = make_stats(counts)

        result = flag_anomalous_overhang_coverage(stats)

        assert result['left_outliers'] == []
        assert result['right_outliers'] == []

    def test_min_contigs_is_configurable(self):
        """Test that the small-assembly floor can be lowered deliberately."""
        counts = {f'contig{i}': (3, 3) for i in range(4)}
        counts['extreme'] = (10_000, 3)
        stats = make_stats(counts)

        result = flag_anomalous_overhang_coverage(stats, min_contigs=5)

        assert result['left_outliers'] == ['extreme']

    def test_outlier_does_not_mask_itself(self):
        """Test that several extreme contigs are all still flagged.

        Mean and standard deviation are both inflated by the outliers being
        looked for, so a cluster of them can hide each other. The median and
        MAD are not affected this way.
        """
        counts = {f'contig{i}': (3, 3) for i in range(12)}
        counts['repeat_a'] = (500, 3)
        counts['repeat_b'] = (480, 3)
        stats = make_stats(counts)

        result = flag_anomalous_overhang_coverage(stats)

        assert set(result['left_outliers']) == {'repeat_a', 'repeat_b'}

    def test_identical_counts_produce_no_flags(self):
        """Test that a zero spread does not divide by zero."""
        stats = make_stats({f'contig{i}': (0, 0) for i in range(10)})

        result = flag_anomalous_overhang_coverage(stats)

        assert result['left_outliers'] == []
        assert result['right_outliers'] == []

    def test_empty_input(self):
        """Test that an empty assembly is handled."""
        result = flag_anomalous_overhang_coverage({})

        assert result['left_outliers'] == []
        assert result['right_outliers'] == []


class TestRankOverhangsByGain:
    """Test the rank_overhangs_by_gain function."""

    def test_rank_overhangs(self):
        """Test ranking overhangs by the sequence they contribute."""
        # All flush with the terminus, so net gain tracks clip length.
        overhangs = [
            OverhangInfo('AT', 2, 1, 100, 'read1', True, 2, 98, 'test_contig', 2),
            OverhangInfo('ATCGATCG', 8, 1, 100, 'read2', True, 8, 92, 'test_contig', 8),
            OverhangInfo('ATCG', 4, 1, 100, 'read3', True, 4, 96, 'test_contig', 4),
        ]

        ranked = rank_overhangs_by_gain(overhangs)

        assert ranked[0].length == 8
        assert ranked[1].length == 4
        assert ranked[2].length == 2

    def test_longest_clip_does_not_win_when_gain_is_smaller(self):
        """Test that a long clip anchored inside the contig loses.

        A 200 bp clip starting 190 bases into the contig contributes only 10
        novel bases, because 190 bases get trimmed before it is grafted on. A
        150 bp clip flush with the terminus contributes all 150. Ranking on raw
        clip length would pick the wrong read.
        """
        buried = OverhangInfo(
            'A' * 200, 200, 191, 1000, 'buried', True, 200, 800, 'test_contig', 10
        )
        flush = OverhangInfo(
            'C' * 150, 150, 1, 1000, 'flush', True, 150, 850, 'test_contig', 150
        )

        ranked = rank_overhangs_by_gain([buried, flush])

        assert ranked[0].read_name == 'flush'
        assert ranked[1].read_name == 'buried'

    def test_least_trim_wins_on_equal_gain(self):
        """Test that the candidate trimming least wins when gains are equal.

        Both contribute 4 novel bases, but one does it without discarding any
        assembly while the other first trims 6 bases.
        """
        intact = OverhangInfo('AAAA', 4, 1, 100, 'intact', True, 4, 96, 'c', 4)
        trimming = OverhangInfo('C' * 10, 10, 7, 100, 'trimming', True, 10, 90, 'c', 4)

        ranked = rank_overhangs_by_gain([intact, trimming])

        assert ranked[0].read_name == 'intact'

    def test_marginally_larger_gain_does_not_justify_trimming(self):
        """Test that one extra base does not buy the right to trim 13.

        Taken from contig01 of the mock dataset. Read L67 contributes one more
        novel base than L91, but only by discarding 13 bases of polished
        consensus and replacing them with a single raw read's sequence. L91
        leaves the contig intact and is the better choice.
        """
        l91 = OverhangInfo('A' * 128, 128, 433, 560, 'L91', False, 128, 122, 'c', 128)
        l67 = OverhangInfo('C' * 142, 142, 440, 547, 'L67', False, 142, 108, 'c', 129)

        ranked = rank_overhangs_by_gain([l67, l91])

        assert ranked[0].read_name == 'L91'

    def test_materially_larger_gain_outranks_an_intact_candidate(self):
        """Test that a gain clearing the margin wins despite trimming.

        Trimming 13 bases to gain 50 more is worth it; trimming 13 to gain 1
        is not. Only the margin separates the two cases.
        """
        intact = OverhangInfo(
            'A' * 128, 128, 433, 560, 'intact', False, 128, 122, 'c', 128
        )
        better = OverhangInfo(
            'C' * 191, 191, 440, 547, 'better', False, 191, 108, 'c', 178
        )

        ranked = rank_overhangs_by_gain([intact, better])

        assert ranked[0].read_name == 'better'

    def test_rank_empty_list(self):
        """Test ranking empty list."""
        ranked = rank_overhangs_by_gain([])
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
        """Test selection of the overhang contributing the most sequence."""
        overhangs = [
            OverhangInfo('AT', 2, 1, 100, 'read1', True, 2, 98, 'test_contig', 2),
            OverhangInfo('ATCGATCG', 8, 1, 100, 'read2', True, 8, 92, 'test_contig', 8),
            OverhangInfo('ATCG', 4, 1, 100, 'read3', True, 4, 96, 'test_contig', 4),
        ]

        best = select_best_overhang(overhangs, min_extension=1, max_homopolymer=50)

        assert best.length == 8
        assert best.read_name == 'read2'

    def test_select_with_min_extension_filter(self):
        """Test selection with minimum extension requirement."""
        overhangs = [
            OverhangInfo('AT', 2, 1, 100, 'read1', True, 2, 98, 'test_contig', 2),
            OverhangInfo('ATCGATCG', 8, 1, 100, 'read2', True, 8, 92, 'test_contig', 8),
            OverhangInfo('ATCG', 4, 1, 100, 'read3', True, 4, 96, 'test_contig', 4),
        ]

        best = select_best_overhang(overhangs, min_extension=5, max_homopolymer=50)

        assert best.length == 8
        assert best.read_name == 'read2'

    def test_min_extension_applies_to_net_gain_not_clip_length(self):
        """Test that min_extension gates on novel sequence contributed.

        A 100 bp clip that only reaches 3 bases past the terminus adds 3 bases,
        so it must not satisfy --min-extension 10 despite its clip length.
        """
        overhangs = [
            OverhangInfo(
                'ATCG' * 25, 100, 98, 1000, 'buried', True, 100, 900, 'test_contig', 3
            ),
        ]

        assert select_best_overhang(overhangs, min_extension=10) is None
        assert select_best_overhang(overhangs, min_extension=3) is not None

    def test_select_with_homopolymer_filtering(self):
        """Test selection avoiding homopolymer runs."""
        overhangs = [
            OverhangInfo(
                'A' * 60, 60, 1, 100, 'read1', True, 60, 40, 'test_contig', 60
            ),  # Has homopolymer
            OverhangInfo(
                'ATCGATCG', 8, 1, 100, 'read2', True, 8, 92, 'test_contig', 8
            ),  # Clean
            OverhangInfo(
                'ATCG', 4, 1, 100, 'read3', True, 4, 96, 'test_contig', 4
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
            OverhangInfo('AT', 2, 1, 100, 'read1', True, 2, 98, 'test_contig', 2),
            OverhangInfo('A', 1, 1, 100, 'read2', True, 1, 99, 'test_contig', 1),
        ]

        best = select_best_overhang(overhangs, min_extension=5, max_homopolymer=50)
        assert best is None


class TestDualEndOverhang:
    """Test cases for reads that overhang both ends of a contig."""

    def setup_method(self):
        """Set up common test data."""
        from teloclip.core.streaming_analysis import process_single_contig_extension

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
            net_gain=18,
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
            net_gain=18,
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

    def test_end_motif_counts_are_windowed_per_end(self):
        """Per-end motif counts cover only the terminal window plus that extension."""
        contig_stats = self.create_dual_end_contig_stats('spanning_read_001')
        # Plant a motif deep in the contig interior; it must not be counted.
        original_sequence = 'A' * 240 + 'TTAGGG' + 'A' * 254

        result = self.process_extension(
            contig_name='short_contig',
            contig_stats=contig_stats,
            original_sequence=original_sequence,
            min_extension=1,
            max_homopolymer=100,
            motif_patterns={'TTAGGG': 'TTAGGG', 'CCCTAA': 'CCCTAA'},
            dry_run=False,
            terminal_length=10,
        )

        assert result is not None

        # Whole-sequence counts still see the interior match.
        assert result.motif_counts['TTAGGG'] == 4  # 3 in the left overhang + 1 interior

        # Per-end windows are 10 + 18 = 28bp and exclude the interior match.
        assert result.end_motif_counts['left'] == {'TTAGGG': 3, 'CCCTAA': 0}
        assert result.end_motif_counts['right'] == {'TTAGGG': 0, 'CCCTAA': 3}

    def test_end_motif_counts_empty_without_patterns(self):
        """No per-end counts are produced when no motifs were requested."""
        contig_stats = self.create_dual_end_contig_stats('spanning_read_001')

        result = self.process_extension(
            contig_name='short_contig',
            contig_stats=contig_stats,
            original_sequence='A' * 500,
            min_extension=1,
            max_homopolymer=100,
            motif_patterns=None,
            dry_run=False,
            terminal_length=10,
        )

        assert result is not None
        assert result.end_motif_counts == {}

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
            net_gain=4,
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
            net_gain=4,
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
            net_gain=60,
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
            net_gain=60,
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
            net_gain=2,
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
            net_gain=2,
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
            net_gain=12,
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
            net_gain=60,
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
            net_gain=24,
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
            net_gain=24,
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
            net_gain=24,
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
            net_gain=24,
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
            net_gain=36,
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
            net_gain=36,
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
                net_gain=8 + 4 * i,
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
                net_gain=8 + 4 * i,
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
