import logging


class _LoggerState:
    """Class to hold the state of the logging infrastructure for the tesserae
    package.

    This class exists to maintain the state of loggers so that we can make them
    formatted nicely with the same state."""

    def __init__(self):
        # Store all our loggers here so we can set levels in them later:
        self._loggers = []

        # Store all our handlers here so we can set name lengths in them later:
        self._handlers = []

        # Keep track of our name length here so we can line up all the logs
        # correctly:
        self._longest_name_length = 12

    @property
    def loggers(self):
        return self._loggers

    @property
    def handlers(self):
        return self._handlers

    @property
    def longest_name_length(self):
        return self._longest_name_length

    @longest_name_length.setter
    def longest_name_length(self, l):
        self._longest_name_length = l

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
_current_logger_state = _LoggerState()


def get_logger(name):
    """Get an all-purpose logger for Tesserae."""

    must_rebuild_formatters = True
    if len(name) > _current_logger_state.longest_name_length:
        _current_logger_state.longest_name_length = len(name)
        must_rebuild_formatters = True

    logger = logging.getLogger(name)
    handler = logging.StreamHandler()
    formatter = _current_logger_state.build_formatter()
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    logger.setLevel(logging.INFO)

    # If we need to rebuild our formatters, we do it here:
    if must_rebuild_formatters:
        _current_logger_state.rebuild_formatters()

    # Add the logger to our list:
    _current_logger_state.add_handler_and_logger(handler, logger)

    return logger


def set_level(lvl):
    """Set the global log level for all loggers we have enabled."""
    _current_logger_state.set_all_levels(lvl)
