import logging
import os
import sys


class CustomFormatter(logging.Formatter):
    """Logging colored formatter, adapted from https://alexandra-zaharia.github.io/posts/make-your-own-custom-color-formatter-with-python-logging"""

    grey = "\x1b[38;21m"
    blue = "\x1b[38;5;39m"
    yellow = "\x1b[38;5;226m"
    red = "\x1b[38;5;196m"
    bold_red = "\x1b[31;1m"
    reset = "\x1b[0m"

    def __init__(self, fmt):
        super().__init__()
        self.fmt = fmt
        self.FORMATS = {
            logging.DEBUG: self.grey + self.fmt + self.reset,
            logging.INFO: self.blue + self.fmt + self.reset,
            logging.WARNING: self.yellow + self.fmt + self.reset,
            logging.ERROR: self.red + self.fmt + self.reset,
            logging.CRITICAL: self.bold_red + self.fmt + self.reset,
        }

    def format(self, record):
        log_fmt = self.FORMATS.get(record.levelno)
        formatter = logging.Formatter(log_fmt)
        return formatter.format(record)


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
