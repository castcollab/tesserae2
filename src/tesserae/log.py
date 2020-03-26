import logging


class _LoggerState:
    """Class to hold the state of the logging infrastructure for the tesserae
    package.

    This class exists to maintain the state of loggers so that we can make them
    formatted nicely with the same state."""

    def __init__(self):
        # Store all our loggers here so we can set levels in them later:
        self.loggers = []

        # Store all our handlers here so we can set name lengths in them later:
        self.handlers = []

        # Keep track of our name length here so we can line up all the logs
        # correctly:
        self.longest_name_length = 12

    def add_handler_and_logger(self, handler, logger):
        self.handlers.append(handler)
        self.loggers.append(logger)

    def build_formatter(self):
        """Build a formatter with space long enough for a name of the given
                length."""
        return logging.Formatter(
            f"%(asctime)s %(name)-{self.longest_name_length}s "
            f"%(levelname)-8s %(message)s"
        )

    def rebuild_formatters(self):
        """Rebuild the formatters for all handlers in `self._handlers`."""
        for handler in self.handlers:
            handler.setFormatter(self.build_formatter())

    def set_all_levels(self, lvl):
        """Set the log level for all loggers."""
        for logger in self.loggers:
            logger.setLevel(lvl)


# Maintain the state of our logs so we can update them all together for nice
# formatting:
_CURRENT_LOGGER_STATE = _LoggerState()


def get_logger(name):
    """Get an all-purpose logger for Tesserae."""

    must_rebuild_formatters = True
    if len(name) > _CURRENT_LOGGER_STATE.longest_name_length:
        _CURRENT_LOGGER_STATE.longest_name_length = len(name)
        must_rebuild_formatters = True

    logger = logging.getLogger(name)
    handler = logging.StreamHandler()
    formatter = _CURRENT_LOGGER_STATE.build_formatter()
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    # If we need to rebuild our formatters, we do it here:
    if must_rebuild_formatters:
        _CURRENT_LOGGER_STATE.rebuild_formatters()

    # Add the logger to our list:
    _CURRENT_LOGGER_STATE.add_handler_and_logger(handler, logger)

    return logger


def set_level(lvl):
    """Set the global log level for all loggers we have enabled."""
    _CURRENT_LOGGER_STATE.set_all_levels(lvl)
