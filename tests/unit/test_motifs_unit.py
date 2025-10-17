"""Unit tests for teloclip.motifs module.

Tests motif pattern generation and sequence matching functions.
"""

from teloclip.motifs import (
    check_sequence_for_patterns,
    construct_regex_pattern,
    count_continuous_runs,
    count_regex_patterns_in_sequence,
    format_pattern_counts,
    make_fuzzy_motif_regex,
    make_motif_regex,
)


class TestMakeMotifRegex:
    """Test exact motif regex generation."""

    def test_make_motif_regex_simple(self):
        """Test simple motif regex generation."""
        result = make_motif_regex('TTAGGG')
        expected = '(TTAGGG)'
        assert result == expected

    def test_make_motif_regex_special_chars(self):
        """Test motif regex with special regex characters."""
        result = make_motif_regex('AT+GC')
        # Special characters should be escaped
        expected = r'(AT\+GC)'
        assert result == expected

    def test_make_motif_regex_empty(self):
        """Test empty motif."""
        result = make_motif_regex('')
        expected = '()'
        assert result == expected


class TestMakeFuzzyMotifRegex:
    """Test fuzzy motif regex generation."""

    def test_make_fuzzy_motif_regex_simple(self):
        """Test fuzzy motif regex for simple repeated patterns."""
        result = make_fuzzy_motif_regex('AAATTT')
        # Should create a pattern that allows for +/- 1 in run lengths
        assert isinstance(result, str)
        assert result.startswith('(')
        assert result.endswith(')')

    def test_make_fuzzy_motif_regex_single_char(self):
        """Test fuzzy motif regex for single characters."""
        result = make_fuzzy_motif_regex('A')
        assert isinstance(result, str)
        assert 'A' in result


class TestCountContinuousRuns:
    """Test continuous run counting."""

    def test_count_continuous_runs_simple(self):
        """Test counting runs in simple sequence."""
        result = count_continuous_runs('AAATTTGGG')
        expected = [('A', 3), ('T', 3), ('G', 3)]
        assert result == expected

    def test_count_continuous_runs_mixed(self):
        """Test counting runs in mixed sequence."""
        result = count_continuous_runs('AATGGCCTA')
        expected = [('A', 2), ('T', 1), ('G', 2), ('C', 2), ('T', 1), ('A', 1)]
        assert result == expected

    def test_count_continuous_runs_single(self):
        """Test counting runs for single character."""
        result = count_continuous_runs('A')
        expected = [('A', 1)]
        assert result == expected

    def test_count_continuous_runs_empty(self):
        """Test counting runs for empty string."""
        result = count_continuous_runs('')
        expected = []
        assert result == expected


class TestConstructRegexPattern:
    """Test regex pattern construction from motif tuples."""

    def test_construct_regex_pattern_simple(self):
        """Test constructing pattern from simple motif tuples."""
        motif_tuples = [('A', 3), ('T', 2)]
        result = construct_regex_pattern(motif_tuples)
        # Should create a pattern allowing +/- 1 in counts
        assert isinstance(result, str)
        assert 'A' in result
        assert 'T' in result

    def test_construct_regex_pattern_single(self):
        """Test constructing pattern from single motif tuple."""
        motif_tuples = [('G', 5)]
        result = construct_regex_pattern(motif_tuples)
        assert isinstance(result, str)
        assert 'G' in result


class TestCountRegexPatternsInSequence:
    """Test counting regex patterns in sequences."""

    def test_count_regex_patterns_simple(self):
        """Test counting simple patterns."""
        sequence = 'TTAGGGTTAGGGTTAGGG'
        patterns = ['TTAGGG', 'CCCTAA']

        result = count_regex_patterns_in_sequence(sequence, patterns)

        assert isinstance(result, dict)
        assert 'TTAGGG' in result
        assert 'CCCTAA' in result
        assert result['TTAGGG'] == 3  # Should find 3 occurrences
        assert result['CCCTAA'] == 0  # Should find 0 occurrences

    def test_count_regex_patterns_overlapping(self):
        """Test counting patterns that might overlap."""
        sequence = 'AAAAAAA'
        patterns = ['AA', 'AAA']

        result = count_regex_patterns_in_sequence(sequence, patterns)

        assert isinstance(result, dict)
        assert 'AA' in result
        assert 'AAA' in result
        # Results depend on regex implementation but should be reasonable
        assert result['AA'] >= 1

    def test_count_regex_patterns_no_matches(self):
        """Test counting patterns with no matches."""
        sequence = 'ATCGATCG'
        patterns = ['TTAGGG', 'CCCTAA']

        result = count_regex_patterns_in_sequence(sequence, patterns)

        assert result['TTAGGG'] == 0
        assert result['CCCTAA'] == 0


class TestFormatPatternCounts:
    """Test pattern count formatting."""

    def test_format_pattern_counts_simple(self):
        """Test formatting simple pattern counts."""
        pattern_counts = {'TTAGGG': 5, 'CCCTAA': 2}
        result = format_pattern_counts(pattern_counts)

        assert isinstance(result, str)
        assert 'TTAGGG' in result
        assert 'CCCTAA' in result
        assert '5' in result
        assert '2' in result

    def test_format_pattern_counts_empty(self):
        """Test formatting empty pattern counts."""
        pattern_counts = {}
        result = format_pattern_counts(pattern_counts)

        assert isinstance(result, str)


class TestCheckSequenceForPatterns:
    """Test sequence pattern checking."""

    def test_check_sequence_for_patterns_match(self):
        """Test checking sequences that should match."""
        sequence = 'TTAGGGTTAGGGTTAGGG'
        patterns = ['TTAGGG']

        result = check_sequence_for_patterns(sequence, patterns)
        assert result is True

    def test_check_sequence_for_patterns_no_match(self):
        """Test checking sequences that should not match."""
        sequence = 'ATCGATCGATCG'
        patterns = ['TTAGGG']

        result = check_sequence_for_patterns(sequence, patterns)
        assert result is False

    def test_check_sequence_for_patterns_multiple_patterns(self):
        """Test checking with multiple patterns."""
        sequence = 'CCCTAACCCTAA'
        patterns = ['TTAGGG', 'CCCTAA']

        result = check_sequence_for_patterns(sequence, patterns)
        assert result is True  # Should match CCCTAA

    def test_check_sequence_for_patterns_min_repeats(self):
        """Test checking with minimum repeat requirements."""
        sequence = 'TTAGGGTTAGGGTTAGGG'  # 3 repeats
        patterns = ['TTAGGG']

        # Should match with low repeat requirement
        result = check_sequence_for_patterns(sequence, patterns, min_repeats=2)
        assert result is True

        # Should not match with high repeat requirement
        result = check_sequence_for_patterns(sequence, patterns, min_repeats=5)
        assert result is False

    def test_check_sequence_for_patterns_edge_cases(self):
        """Test edge cases for pattern checking."""
        # Empty sequence
        result = check_sequence_for_patterns('', ['TTAGGG'])
        assert result is False

        # Empty patterns list
        result = check_sequence_for_patterns('TTAGGG', [])
        assert result is False

        # Single character pattern
        result = check_sequence_for_patterns('AAAA', ['A'])
        assert result is True


class TestMotifIntegration:
    """Test integration between different motif functions."""

    def test_motif_workflow(self):
        """Test a complete motif analysis workflow."""
        # Create some motif patterns
        exact_pattern = make_motif_regex('TTAGGG')
        fuzzy_pattern = make_fuzzy_motif_regex('AAATTT')

        # Test sequence
        sequence = 'TTAGGGTTAGGGAAATTTAAATTT'

        # Check if patterns match
        exact_match = check_sequence_for_patterns(sequence, [exact_pattern])
        fuzzy_match = check_sequence_for_patterns(sequence, [fuzzy_pattern])

        assert exact_match is True
        # Fuzzy match depends on implementation but should not crash
        assert isinstance(fuzzy_match, bool)

    def test_pattern_counting_workflow(self):
        """Test pattern counting workflow."""
        # Generate patterns
        patterns = [make_motif_regex('TTAGGG'), make_motif_regex('CCCTAA')]

        # Test sequence
        sequence = 'TTAGGGTTAGGGCCCTAACCCTAA'

        # Count patterns
        counts = count_regex_patterns_in_sequence(sequence, patterns)

        # Format results
        formatted = format_pattern_counts(counts)

        assert isinstance(counts, dict)
        assert isinstance(formatted, str)
        assert len(counts) == len(patterns)
