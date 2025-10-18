"""
Logging configuration for the teloclip package.

This module provides functionality to initialize and configure logging with rich
formatting for better readability in terminal output. It uses the 'rich' library
to create visually enhanced log messages.
"""

import logging
from pathlib import Path
from typing import Optional, Union

from rich.console import Console
from rich.logging import RichHandler


def init_logging(
    loglevel: str = 'DEBUG', logfile: Optional[Union[str, Path]] = None
) -> None:
    """
    Initialize root logger with specified log level and rich formatting.

    Configures the global logging system with rich formatting for console output
    and optionally writes logs to a file.

    Parameters
    ----------
    loglevel : str, optional
        The log level to use (e.g., "DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"),
        by default "DEBUG".
    logfile : str or Path, optional
        Path to log file. If provided, logs will be written to both console and file.

    Returns
    -------
    None
        This function configures the global logging system and doesn't return a value.

    Raises
    ------
    ValueError
        If the provided log level is invalid.
    OSError
        If the log file cannot be created or written to.
    """
    # Convert log level string to numeric value
    numeric_level = getattr(logging, loglevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError(f'Invalid log level: {loglevel}')

    # Get the root logger
    root_logger = logging.getLogger()

    # Clear existing handlers if any are present
    if root_logger.hasHandlers():
        root_logger.handlers.clear()

    # Set the logger's level according to specified loglevel
    root_logger.setLevel(numeric_level)

    # Add rich handler for console output (stderr for better visibility in scripts)
    console_handler = RichHandler(console=Console(stderr=True))
    console_handler.setLevel(numeric_level)
    root_logger.addHandler(console_handler)

    # Add file handler if logfile is specified
    if logfile:
        try:
            logfile_path = Path(logfile)

            # Create parent directories if they don't exist
            logfile_path.parent.mkdir(parents=True, exist_ok=True)

            # Create file handler with detailed formatting
            file_handler = logging.FileHandler(logfile_path, mode='w', encoding='utf-8')
            file_handler.setLevel(numeric_level)

            # Create detailed formatter for file output
            file_formatter = logging.Formatter(
                fmt='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                datefmt='%Y-%m-%d %H:%M:%S',
            )
            file_handler.setFormatter(file_formatter)

            root_logger.addHandler(file_handler)

            # Log the initialization
            logging.info(f'Logging initialized with level: {loglevel}')
            logging.info(f'Log file: {logfile_path.absolute()}')

        except OSError as e:
            # If file logging fails, continue with console-only logging
            logging.warning(f'Failed to create log file {logfile}: {e}')
            logging.warning('Continuing with console-only logging')
    else:
        logging.info(f'Logging initialized with level: {loglevel}')
