"""Unit tests for teloclip.commands.extend parsing functions.

Tests the contig exclusion and motif parsing utilities.
"""

from pathlib import Path
import re
import tempfile
from unittest.mock import patch

import click
import pytest

from teloclip.commands.extend import (
    combine_excluded_contigs,
    get_motif_regex,
    parse_excluded_contigs,
    read_excluded_contigs_file,
)


class TestReadExcludedContigsFile:
    """Test reading contig names from exclusion files."""

    def test_read_excluded_contigs_file_basic(self):
        """Test reading a basic exclusion file."""
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as f:
            f.write('contig1\ncontig2\ncontig3\n')
            f.flush()
            temp_path = Path(f.name)

        try:
            contigs = read_excluded_contigs_file(temp_path)
            assert contigs == ['contig1', 'contig2', 'contig3']
        finally:
            temp_path.unlink()

    def test_read_excluded_contigs_file_with_whitespace(self):
        """Test reading file with various whitespace."""
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as f:
            f.write('  contig1  \n\tcontig2\t\n  \n contig3 \n\n')
            f.flush()
            temp_path = Path(f.name)

        try:
            contigs = read_excluded_contigs_file(temp_path)
            assert contigs == ['contig1', 'contig2', 'contig3']
        finally:
            temp_path.unlink()

    def test_read_excluded_contigs_file_different_line_endings(self):
        """Test reading file with different line ending styles."""
        # Test with Windows-style line endings
        with tempfile.NamedTemporaryFile(mode='wb', delete=False, suffix='.txt') as f:
            f.write(b'contig1\r\ncontig2\r\ncontig3\r\n')
            temp_path = Path(f.name)

        try:
            contigs = read_excluded_contigs_file(temp_path)
            assert contigs == ['contig1', 'contig2', 'contig3']
        finally:
            temp_path.unlink()

    def test_read_excluded_contigs_file_empty_file(self):
        """Test error handling for empty file."""
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as f:
            f.write('')
            temp_path = Path(f.name)

        try:
            with pytest.raises(click.ClickException) as exc_info:
                read_excluded_contigs_file(temp_path)
            assert 'No valid contig names found' in str(exc_info.value)
        finally:
            temp_path.unlink()

    def test_read_excluded_contigs_file_only_empty_lines(self):
        """Test error handling for file with only empty lines."""
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as f:
            f.write('\n\n  \n\t\n   \n')
            temp_path = Path(f.name)

        try:
            with pytest.raises(click.ClickException) as exc_info:
                read_excluded_contigs_file(temp_path)
            assert 'No valid contig names found' in str(exc_info.value)
        finally:
            temp_path.unlink()

    def test_read_excluded_contigs_file_nonexistent(self):
        """Test error handling for nonexistent file."""
        nonexistent_path = Path('/nonexistent/path/file.txt')

        with pytest.raises(click.ClickException) as exc_info:
            read_excluded_contigs_file(nonexistent_path)
        assert 'Exclusion file not found' in str(exc_info.value)

    def test_read_excluded_contigs_file_permission_error(self):
        """Test error handling for permission errors."""
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as f:
            f.write('contig1\n')
            temp_path = Path(f.name)

        try:
            # Mock permission error
            with patch(
                'builtins.open', side_effect=PermissionError('Permission denied')
            ):
                with pytest.raises(click.ClickException) as exc_info:
                    read_excluded_contigs_file(temp_path)
                assert 'Error reading exclusion file' in str(exc_info.value)
        finally:
            temp_path.unlink()


class TestParseExcludedContigs:
    """Test parsing excluded contig names from command line string."""

    def test_parse_excluded_contigs_basic(self):
        """Test basic contig parsing."""
        contig_dict = {'contig1': 1000, 'contig2': 2000, 'contig3': 3000}
        exclude_str = 'contig1,contig2'

        result = parse_excluded_contigs(exclude_str, contig_dict)
        assert result == {'contig1', 'contig2'}

    def test_parse_excluded_contigs_with_whitespace(self):
        """Test parsing with extra whitespace."""
        contig_dict = {'contig1': 1000, 'contig2': 2000, 'contig3': 3000}
        exclude_str = ' contig1 , contig2 ,  contig3  '

        result = parse_excluded_contigs(exclude_str, contig_dict)
        assert result == {'contig1', 'contig2', 'contig3'}

    def test_parse_excluded_contigs_empty_string(self):
        """Test parsing empty string."""
        contig_dict = {'contig1': 1000}
        exclude_str = ''

        result = parse_excluded_contigs(exclude_str, contig_dict)
        assert result == set()

    def test_parse_excluded_contigs_none_string(self):
        """Test parsing None string."""
        contig_dict = {'contig1': 1000}
        exclude_str = None

        result = parse_excluded_contigs(exclude_str, contig_dict)
        assert result == set()

    def test_parse_excluded_contigs_invalid_names(self):
        """Test parsing with invalid contig names."""
        contig_dict = {'contig1': 1000, 'contig2': 2000}
        exclude_str = 'contig1,invalid_contig,contig2'

        with patch('logging.warning') as mock_warning:
            result = parse_excluded_contigs(exclude_str, contig_dict)
            assert result == {'contig1', 'contig2'}
            mock_warning.assert_called()

    def test_parse_excluded_contigs_all_invalid(self):
        """Test parsing with all invalid contig names."""
        contig_dict = {'contig1': 1000, 'contig2': 2000}
        exclude_str = 'invalid1,invalid2'

        with patch('logging.warning') as mock_warning:
            result = parse_excluded_contigs(exclude_str, contig_dict)
            assert result == set()
            mock_warning.assert_called()

    def test_parse_excluded_contigs_empty_items(self):
        """Test parsing with empty items in comma-separated list."""
        contig_dict = {'contig1': 1000, 'contig2': 2000}
        exclude_str = 'contig1,,contig2,,'

        result = parse_excluded_contigs(exclude_str, contig_dict)
        assert result == {'contig1', 'contig2'}


class TestCombineExcludedContigs:
    """Test combining excluded contigs from multiple sources."""

    def setup_method(self):
        """Set up test contig dictionary."""
        self.contig_dict = {
            'contig1': 1000,
            'contig2': 2000,
            'contig3': 3000,
            'contig4': 4000,
        }

    def test_combine_excluded_contigs_string_only(self):
        """Test combining with only string source."""
        exclude_str = 'contig1,contig2'
        exclude_file = None

        result = combine_excluded_contigs(exclude_str, exclude_file, self.contig_dict)
        assert result == {'contig1', 'contig2'}

    def test_combine_excluded_contigs_file_only(self):
        """Test combining with only file source."""
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as f:
            f.write('contig3\ncontig4\n')
            temp_path = Path(f.name)

        try:
            exclude_str = None
            result = combine_excluded_contigs(exclude_str, temp_path, self.contig_dict)
            assert result == {'contig3', 'contig4'}
        finally:
            temp_path.unlink()

    def test_combine_excluded_contigs_both_sources(self):
        """Test combining from both string and file sources."""
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as f:
            f.write('contig3\ncontig4\n')
            temp_path = Path(f.name)

        try:
            exclude_str = 'contig1,contig2'
            with patch('logging.warning') as mock_warning:
                result = combine_excluded_contigs(
                    exclude_str, temp_path, self.contig_dict
                )
                assert result == {'contig1', 'contig2', 'contig3', 'contig4'}
                # Should warn about using both sources
                mock_warning.assert_called()
        finally:
            temp_path.unlink()

    def test_combine_excluded_contigs_duplicates(self):
        """Test handling of duplicate contig names."""
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as f:
            f.write('contig1\ncontig3\n')
            temp_path = Path(f.name)

        try:
            exclude_str = 'contig1,contig2'  # contig1 appears in both
            with patch('logging.info') as mock_info:
                result = combine_excluded_contigs(
                    exclude_str, temp_path, self.contig_dict
                )
                assert result == {'contig1', 'contig2', 'contig3'}
                # Should log about removing duplicates
                mock_info.assert_called()
        finally:
            temp_path.unlink()

    def test_combine_excluded_contigs_no_sources(self):
        """Test combining with no sources provided."""
        exclude_str = None
        exclude_file = None

        result = combine_excluded_contigs(exclude_str, exclude_file, self.contig_dict)
        assert result == set()

    def test_combine_excluded_contigs_invalid_names(self):
        """Test combining with invalid contig names."""
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as f:
            f.write('invalid_contig\ncontig2\n')
            temp_path = Path(f.name)

        try:
            exclude_str = 'contig1,another_invalid'
            with patch('logging.warning') as mock_warning:
                result = combine_excluded_contigs(
                    exclude_str, temp_path, self.contig_dict
                )
                assert result == {'contig1', 'contig2'}
                mock_warning.assert_called()
        finally:
            temp_path.unlink()


class TestGetMotifRegex:
    """Test motif regex pattern generation."""

    def test_get_motif_regex_basic(self):
        """Test basic motif regex generation."""
        motif_str = 'TTAGGG'

        result = get_motif_regex(motif_str, fuzzy=False)

        # Should contain both forward and reverse complement
        assert 'TTAGGG' in result
        assert 'CCCTAA' in result  # reverse complement
        assert len(result) == 2

        # Check patterns are compiled regex objects
        assert isinstance(result['TTAGGG'], re.Pattern)
        assert isinstance(result['CCCTAA'], re.Pattern)

    def test_get_motif_regex_multiple_motifs(self):
        """Test regex generation for multiple motifs."""
        motif_str = 'TTAGGG,TTAATTGG'

        result = get_motif_regex(motif_str, fuzzy=False)

        # Should contain all motifs and their reverse complements
        expected_motifs = {'TTAGGG', 'CCCTAA', 'TTAATTGG', 'CCAATTAA'}
        assert set(result.keys()) == expected_motifs

    def test_get_motif_regex_with_whitespace(self):
        """Test motif parsing with extra whitespace."""
        motif_str = ' TTAGGG , TTAATTGG '

        result = get_motif_regex(motif_str, fuzzy=False)

        expected_motifs = {'TTAGGG', 'CCCTAA', 'TTAATTGG', 'CCAATTAA'}
        assert set(result.keys()) == expected_motifs

    def test_get_motif_regex_fuzzy_mode(self):
        """Test fuzzy motif regex generation."""
        motif_str = 'TTAGGG'

        result = get_motif_regex(motif_str, fuzzy=True)

        # In fuzzy mode, pattern names should have '_fuzzy' suffix
        assert 'TTAGGG_fuzzy' in result
        assert 'CCCTAA_fuzzy' in result
        assert len(result) == 2

    def test_get_motif_regex_invalid_bases(self):
        """Test handling of invalid DNA bases."""
        motif_str = 'TTAGGG,INVALID123,CCCTAA'

        with patch('logging.warning') as mock_warning:
            result = get_motif_regex(motif_str, fuzzy=False)

            # Should only contain valid motifs
            valid_keys = [k for k in result.keys() if not k.endswith('_fuzzy')]
            expected_valid = {'TTAGGG', 'CCCTAA'}  # INVALID123 should be filtered out
            assert set(valid_keys) == expected_valid

            # Should log warning about invalid motif
            mock_warning.assert_called()

    def test_get_motif_regex_empty_motifs(self):
        """Test handling of empty motifs in list."""
        motif_str = 'TTAGGG,,CCCTAA,'

        result = get_motif_regex(motif_str, fuzzy=False)

        # Should ignore empty motifs
        expected_motifs = {'TTAGGG', 'CCCTAA'}
        actual_motifs = set(result.keys())
        assert actual_motifs == expected_motifs

    def test_get_motif_regex_all_invalid(self):
        """Test handling when all motifs are invalid."""
        motif_str = 'INVALID123,XYZ789'

        with patch('logging.warning') as mock_warning:
            result = get_motif_regex(motif_str, fuzzy=False)

            # Should return empty dictionary
            assert result == {}
            mock_warning.assert_called()

    def test_get_motif_regex_case_insensitive(self):
        """Test that motifs are converted to uppercase."""
        motif_str = 'ttaggg,ccctaa'

        result = get_motif_regex(motif_str, fuzzy=False)

        # Should be converted to uppercase
        expected_motifs = {'TTAGGG', 'CCCTAA'}
        assert set(result.keys()) == expected_motifs

    def test_get_motif_regex_duplicate_motifs(self):
        """Test handling of duplicate motifs."""
        motif_str = 'TTAGGG,TTAGGG,CCCTAA'

        result = get_motif_regex(motif_str, fuzzy=False)

        # Should deduplicate - TTAGGG and CCCTAA are reverse complements
        expected_motifs = {'TTAGGG', 'CCCTAA'}
        assert set(result.keys()) == expected_motifs

    def test_get_motif_regex_self_complement(self):
        """Test motif that is its own reverse complement."""
        motif_str = 'ATAT'  # ATAT reverse complement is ATAT

        result = get_motif_regex(motif_str, fuzzy=False)

        # Should only have one entry since it's self-complementary
        assert set(result.keys()) == {'ATAT'}
        assert len(result) == 1
