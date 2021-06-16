"""Handle coloring for logger

Note: all non-Windows platforms are supporting ANSI escapes so we use them
"""

import logging


def log_coloring(fn):
    """Handler for log coloring

    Based on https://stackoverflow.com/questions/384076/how-can-i-color-python-logging-output
    Color code explanation here: http://jafrog.com/2013/11/23/colors-in-terminal.html
    """

    def new(*args):
        levelno = args[1].levelno
        if levelno >= logging.CRITICAL:
            color = "\x1b[31m"  # red
        elif levelno >= logging.ERROR:
            color = "\x1b[31m"  # red
        elif levelno >= logging.WARNING:
            color = "\x1b[93m"  # light green
        elif levelno >= logging.INFO:
            color = "\x1b[32m"  # light yellow
        elif levelno >= logging.DEBUG:
            color = "\x1b[36m"  # pink
        else:
            color = "\x1b[0m"  # normal
        args[1].msg = color + args[1].msg + "\x1b[0m"
        return fn(*args)

    return new


logging.StreamHandler.emit = log_coloring(logging.StreamHandler.emit)
logger = logging.getLogger("nmdb")
formatter = logging.Formatter(
    "[%(asctime)s] [%(levelname)8s] --- %(message)-50s (%(filename)s:%(lineno)s)",
    "%Y-%m-%d %H:%M:%S",
)
stream = logging.StreamHandler()
stream.setFormatter(formatter)
logger.setLevel(logging.INFO)
logger.addHandler(stream)
