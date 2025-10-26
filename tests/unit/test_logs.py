"""Unit tests for teloclip.logs module.

Tests logging initialization and configuration functions.
"""

import logging
from pathlib import Path
import tempfile
from unittest.mock import patch

import pytest
from rich.logging import RichHandler

from teloclip.logs import init_logging


class TestInitLogging:
    """Test logging initialization function."""

    def teardown_method(self):
        """Clean up logging configuration after each test."""
        # Clear all handlers from root logger
        root_logger = logging.getLogger()
        root_logger.handlers.clear()
        root_logger.setLevel(logging.WARNING)  # Reset to default

    def test_init_logging_debug_level(self):
        """Test logging initialization with DEBUG level."""
        init_logging('DEBUG')

        root_logger = logging.getLogger()
        assert root_logger.level == logging.DEBUG
        assert len(root_logger.handlers) == 1
        assert isinstance(root_logger.handlers[0], RichHandler)

    def test_init_logging_info_level(self):
        """Test logging initialization with INFO level."""
        init_logging('INFO')

        root_logger = logging.getLogger()
        assert root_logger.level == logging.INFO
        assert len(root_logger.handlers) == 1
        assert isinstance(root_logger.handlers[0], RichHandler)

    def test_init_logging_warning_level(self):
        """Test logging initialization with WARNING level."""
        init_logging('WARNING')

        root_logger = logging.getLogger()
        assert root_logger.level == logging.WARNING
        assert len(root_logger.handlers) == 1

    def test_init_logging_error_level(self):
        """Test logging initialization with ERROR level."""
        init_logging('ERROR')

        root_logger = logging.getLogger()
        assert root_logger.level == logging.ERROR

    def test_init_logging_critical_level(self):
        """Test logging initialization with CRITICAL level."""
        init_logging('CRITICAL')

        root_logger = logging.getLogger()
        assert root_logger.level == logging.CRITICAL

    def test_init_logging_case_insensitive(self):
        """Test that log level is case insensitive."""
        init_logging('info')
        assert logging.getLogger().level == logging.INFO

        init_logging('Warning')
        assert logging.getLogger().level == logging.WARNING

    def test_init_logging_invalid_level(self):
        """Test error handling for invalid log level."""
        with pytest.raises(ValueError) as exc_info:
            init_logging('INVALID_LEVEL')

        assert 'Invalid log level: INVALID_LEVEL' in str(exc_info.value)

    def test_init_logging_with_logfile(self):
        """Test logging initialization with file output."""
        with tempfile.TemporaryDirectory() as temp_dir:
            logfile = Path(temp_dir) / 'test.log'

            init_logging('INFO', logfile)

            root_logger = logging.getLogger()
            assert len(root_logger.handlers) == 2  # Console + file

            # Check that file handler is added
            file_handlers = [
                h for h in root_logger.handlers if isinstance(h, logging.FileHandler)
            ]
            assert len(file_handlers) == 1

            # Check file was created
            assert logfile.exists()

    def test_init_logging_with_nested_logfile_path(self):
        """Test logging with logfile in nested directory."""
        with tempfile.TemporaryDirectory() as temp_dir:
            logfile = Path(temp_dir) / 'subdir' / 'nested' / 'test.log'

            init_logging('DEBUG', logfile)

            # Check that nested directories were created
            assert logfile.parent.exists()
            assert logfile.exists()

    def test_init_logging_file_permission_error(self):
        """Test handling of file permission errors."""
        # Use a path that should cause permission issues
        invalid_path = '/root/cannot_write_here.log'

        with patch('logging.FileHandler') as mock_file_handler:
            mock_file_handler.side_effect = PermissionError('Permission denied')

            # Should not raise exception, should continue with console logging
            init_logging('INFO', invalid_path)

            root_logger = logging.getLogger()
            # Should only have console handler due to file error
            assert len(root_logger.handlers) == 1
            assert isinstance(root_logger.handlers[0], RichHandler)

    def test_init_logging_clears_existing_handlers(self):
        """Test that existing handlers are cleared."""
        # Add a dummy handler first
        root_logger = logging.getLogger()
        dummy_handler = logging.StreamHandler()
        root_logger.addHandler(dummy_handler)

        init_logging('INFO')

        # Should only have the new RichHandler
        assert len(root_logger.handlers) == 1
        assert isinstance(root_logger.handlers[0], RichHandler)
        assert dummy_handler not in root_logger.handlers

    def test_init_logging_default_level(self):
        """Test logging initialization with default DEBUG level."""
        init_logging()

        root_logger = logging.getLogger()
        assert root_logger.level == logging.DEBUG

    def test_init_logging_file_formatter(self):
        """Test that file handler gets proper formatter."""
        with tempfile.TemporaryDirectory() as temp_dir:
            logfile = Path(temp_dir) / 'test.log'

            init_logging('INFO', logfile)

            root_logger = logging.getLogger()
            file_handlers = [
                h for h in root_logger.handlers if isinstance(h, logging.FileHandler)
            ]

            assert len(file_handlers) == 1
            file_handler = file_handlers[0]
            assert file_handler.formatter is not None

            # Check formatter format string contains expected elements
            fmt = file_handler.formatter._fmt
            assert '%(asctime)s' in fmt
            assert '%(name)s' in fmt
            assert '%(levelname)s' in fmt
            assert '%(message)s' in fmt

    def test_init_logging_with_pathlib_path(self):
        """Test logging initialization with Path object."""
        with tempfile.TemporaryDirectory() as temp_dir:
            logfile = Path(temp_dir) / 'test.log'

            init_logging('INFO', logfile)

            root_logger = logging.getLogger()
            assert len(root_logger.handlers) == 2
            assert logfile.exists()
