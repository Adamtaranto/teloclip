"""Unit tests for teloclip.utils module.

Tests utility functions.
"""

from unittest.mock import patch

import pytest

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
