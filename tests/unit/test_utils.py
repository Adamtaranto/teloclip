"""Unit tests for teloclip.utils module.

Tests utility functions.
"""

from pathlib import Path
import sys
import tempfile
from unittest.mock import Mock, patch

import click
import pytest

from teloclip.commands.extend import validate_output_directories
from teloclip.extract_io import create_fasta_index
from teloclip.utils import isfile


class TestIsFile:
    """Test file existence checking function."""

    @patch('os.path.abspath')
    @patch('os.path.isfile')
    def test_isfile_exists(self, mock_isfile, mock_abspath):
        """Test checking for existing file."""
        mock_isfile.return_value = True
        mock_abspath.return_value = '/absolute/path/to/existing_file.txt'

        result = isfile('existing_file.txt')
        assert result == '/absolute/path/to/existing_file.txt'
        mock_isfile.assert_called_once_with('existing_file.txt')
        mock_abspath.assert_called_once_with('existing_file.txt')

    @patch('os.path.isfile')
    def test_isfile_not_exists(self, mock_isfile):
        """Test checking for non-existing file (should exit)."""
        mock_isfile.return_value = False

        with pytest.raises(SystemExit):
            isfile('non_existing_file.txt')
        mock_isfile.assert_called_once_with('non_existing_file.txt')

    @patch('os.path.isfile')
    def test_isfile_empty_path(self, mock_isfile):
        """Test checking empty path (should exit)."""
        mock_isfile.return_value = False

        with pytest.raises(SystemExit):
            isfile('')
        mock_isfile.assert_called_once_with('')


class TestValidateOutputDirectories:
    """Test output directory validation and creation."""

    def test_validate_output_directories_existing_dirs(self):
        """Test validation when directories already exist."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            output_fasta = temp_path / 'existing' / 'output.fasta'
            stats_report = temp_path / 'existing2' / 'report.txt'

            # Create the directories
            output_fasta.parent.mkdir(parents=True)
            stats_report.parent.mkdir(parents=True)

            # Should not raise exception
            validate_output_directories(output_fasta, stats_report)

    def test_validate_output_directories_create_needed_dirs(self):
        """Test creating directories when they don't exist."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            output_fasta = temp_path / 'new_dir1' / 'subdir' / 'output.fasta'
            stats_report = temp_path / 'new_dir2' / 'subdir' / 'report.txt'

            # Directories don't exist yet
            assert not output_fasta.parent.exists()
            assert not stats_report.parent.exists()

            validate_output_directories(output_fasta, stats_report)

            # Should create the directories
            assert output_fasta.parent.exists()
            assert stats_report.parent.exists()

    def test_validate_output_directories_none_paths(self):
        """Test handling of None paths."""
        # Should not raise exception with None paths
        validate_output_directories(None, None)

    def test_validate_output_directories_mixed_none(self):
        """Test handling with one None path."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            output_fasta = temp_path / 'new_dir' / 'output.fasta'

            validate_output_directories(output_fasta, None)

            # Should create directory for non-None path
            assert output_fasta.parent.exists()

    def test_validate_output_directories_permission_error(self):
        """Test handling of permission errors during directory creation."""
        invalid_path = Path('/some/path/output.fasta')
        with patch.object(
            Path, 'mkdir', side_effect=PermissionError('Permission denied')
        ):
            with pytest.raises(click.ClickException) as exc_info:
                validate_output_directories(invalid_path, None)
        assert 'Cannot create output directory' in str(exc_info.value)


class TestCreateFastaIndex:
    """Test FASTA index creation."""

    def test_create_fasta_index_already_exists(self):
        """Test when index already exists."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            fasta_file = temp_path / 'test.fasta'
            index_file = temp_path / 'test.fasta.fai'

            # Create both files
            fasta_file.touch()
            index_file.touch()

            # Mock pysam module that gets imported within the function
            with patch.dict('sys.modules', {'pysam': Mock()}):
                with patch('pysam.faidx') as mock_faidx:
                    result = create_fasta_index(fasta_file)

                    # Should return existing index path without calling pysam.faidx
                    assert result == index_file
                    mock_faidx.assert_not_called()

    def test_create_fasta_index_needs_creation(self):
        """Test creating index when it doesn't exist."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            fasta_file = temp_path / 'test.fasta'
            index_file = temp_path / 'test.fasta.fai'

            # Create only FASTA file
            fasta_file.touch()
            assert not index_file.exists()

            with patch.dict('sys.modules', {'pysam': Mock()}):
                with patch('pysam.faidx') as mock_faidx:
                    with patch('builtins.print'):  # Suppress print output
                        result = create_fasta_index(fasta_file)

                    # Should call pysam.faidx and return index path
                    assert result == index_file
                    mock_faidx.assert_called_once_with(str(fasta_file))

    def test_create_fasta_index_prints_message(self):
        """Test that appropriate message is printed when creating index."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            fasta_file = temp_path / 'test.fasta'
            index_file = temp_path / 'test.fasta.fai'

            fasta_file.touch()

            with patch.dict('sys.modules', {'pysam': Mock()}):
                with patch('pysam.faidx'):
                    with patch('builtins.print') as mock_print:
                        create_fasta_index(fasta_file)

                        # Should print message about creating index
                        mock_print.assert_called_once()
                        args = mock_print.call_args[0]
                        assert 'Creating FASTA index' in args[0]
                        assert str(index_file) in args[0]
                        # Check that it prints to stderr
                        assert mock_print.call_args[1]['file'] == sys.stderr
