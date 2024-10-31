import logging
import os
import sys


def isfile(path):
    """
    Check if file exists.

    Args:
        path (_str_): Path to file

    Returns:
        _str_: Returns abspath to file if found, else prints warning and exits.
    """
    if not os.path.isfile(path):
        logging.error("Input file not found: %s" % path)
        sys.exit(1)
    else:
        return os.path.abspath(path)
