"""Unit tests for teloclip.utils module.

Tests utility functions.
"""

from unittest.mock import patch
from teloclip.utils import isfile


class TestIsFile:
    """Test file existence checking function."""

    @patch('os.path.isfile')
    def test_isfile_exists(self, mock_isfile):
        """Test checking for existing file."""
        mock_isfile.return_value = True

        result = isfile('existing_file.txt')
        assert result is True
        mock_isfile.assert_called_once_with('existing_file.txt')

    @patch('os.path.isfile')
    def test_isfile_not_exists(self, mock_isfile):
        """Test checking for non-existing file."""
        mock_isfile.return_value = False

        result = isfile('non_existing_file.txt')
        assert result is False
        mock_isfile.assert_called_once_with('non_existing_file.txt')

    @patch('os.path.isfile')
    def test_isfile_empty_path(self, mock_isfile):
        """Test checking empty path."""
        mock_isfile.return_value = False

        result = isfile('')
        assert result is False
        mock_isfile.assert_called_once_with('')
