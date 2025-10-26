"""
Utility functions for teloclip.

This module provides common utility functions used throughout the teloclip
package, including file validation and path handling utilities.
"""

import logging
import os
import sys


def isfile(path):
    """
    Check if file exists and return absolute path.

    Parameters
    ----------
    path : str
        Path to file to check.

    Returns
    -------
    str
        Absolute path to file if found.

    Raises
    ------
    SystemExit
        If file is not found, logs error and exits with code 1.
    """
    if not os.path.isfile(path):
        logging.error('Input file not found: %s' % path)
        sys.exit(1)
    else:
        return os.path.abspath(path)
