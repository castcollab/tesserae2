import logging

# Store all our loggers here so we can set levels in them later:
__loggers = []

# Store all our handlers here so we can set name lengths in them later:
__handlers = []

# Redefine log levels here so we can have more granular control later (if we want):
CRITICAL = logging.CRITICAL
ERROR = logging.ERROR
WARNING = logging.WARNING
INFO = logging.INFO
DEBUG = logging.DEBUG
NOTSET = logging.NOTSET

# Keep track of our name length here so we can line up all the logs correctly:
__longest_name_length = 12


def __add_handler_and_logger(handler, logger):
    """Add the given handler and logger to our internal lists so we can keep track of them."""
    __handlers.append(handler)
    __loggers.append(logger)


def _build_formatter(name_length):
    """Build a formatter with space long enough for a name of the given length."""
    return logging.Formatter(
        "%(asctime)s %(name)-" + str(name_length) + "s %(levelname)-8s %(message)s"
    )


def get_logger(name):
    """Get an all-purpose logger for Tesserae."""

    # Not sure of a good way around this:
    global __loggers, __handlers, __longest_name_length

    must_rebuild_formatters = True
    if len(name) > __longest_name_length:
        __longest_name_length = len(name) + 2
        must_rebuild_formatters = True

    logger = logging.getLogger(name)
    handler = logging.StreamHandler()
    formatter = _build_formatter(__longest_name_length)
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    logger.setLevel(logging.INFO)

    # If we need to rebuild our formatters, we do it here:
    if must_rebuild_formatters:
        for handler in __handlers:
            handler.setFormatter(_build_formatter(__longest_name_length))

    # Add the logger to our list:
    __add_handler_and_logger(handler, logger)

    return logger


def set_level(lvl):
    """Set the global log level for all loggers we have enabled."""
    for logger in __loggers:
        logger.setLevel(lvl)
