import logging
import sys


def init_logging():
    # Define the log message format
    fmt = "%(asctime)s | %(levelname)s | %(module)s:%(lineno)s:%(funcName)s() | %(message)s"

    # Set up a StreamHandler for stderr instead of stdout
    handler_sh = logging.StreamHandler(sys.stderr)
    handler_sh.setFormatter(CustomFormatter(fmt))

    # Configure the logging with our handler for stderr
    logging.basicConfig(format=fmt, level=logging.DEBUG, handlers=[handler_sh])


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
