"""Unit tests for teloclip.core.motifs module.

Tests motif pattern generation and sequence matching functions.
"""

from teloclip.core.motifs import (
    check_sequence_for_patterns,
    construct_regex_pattern,
    count_continuous_runs,
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
        # Each run becomes a +/- 1 quantified range.
        assert result == r'A{2,4}T{1,3}'

    def test_construct_regex_pattern_single(self):
        """Test constructing pattern from single motif tuple."""
        motif_tuples = [('G', 5)]
        result = construct_regex_pattern(motif_tuples)
        assert result == r'G{4,6}'

    def test_construct_regex_pattern_telomeric_motif(self):
        """A run of one stays literal, so only the real homopolymers flex."""
        motif_tuples = count_continuous_runs('TTTAGGG')
        assert construct_regex_pattern(motif_tuples) == r'T{2,4}AG{2,4}'

    def test_construct_regex_pattern_no_repeated_bases(self):
        """A motif with no runs longer than one is reproduced verbatim."""
        motif_tuples = count_continuous_runs('ATATAT')
        assert construct_regex_pattern(motif_tuples) == r'ATATAT'

    def test_construct_regex_pattern_multiple_runs(self):
        """Every run is quantified independently."""
        motif_tuples = count_continuous_runs('AAAAGGGTTTCCCC')
        assert construct_regex_pattern(motif_tuples) == r'A{3,5}G{2,4}T{2,4}C{3,5}'

    def test_construct_regex_pattern_empty(self):
        """An empty motif produces an empty pattern."""
        assert construct_regex_pattern([]) == ''


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
