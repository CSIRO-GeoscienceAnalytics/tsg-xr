
import logging
from pathlib import Path


def rm_tree(pth: Path):
    """Recursively remove a folder."""
    if pth.exists():
        if pth.is_file():
            pth.unlink()
        else:
            for child in pth.iterdir():
                if child.is_file():
                    child.unlink()
                else:
                    rm_tree(child)
            pth.rmdir()


def Handle(
    logger,
    handler_class=logging.StreamHandler,
    formatter="%(asctime)s %(name)s - %(levelname)s: %(message)s",
    level=None,
):
    """
    Handle a logger with a standardised formatting.

    Parameters
    -----------
    logger : :class:`logging.Logger` | :class:`str`
        Logger or module name to source a logger from.
    handler_class : :class:`logging.Handler`
        Handler class for the logging messages.
    formatter : :class:`str` | :class:`logging.Formatter`
        Formatter for the logging handler. Strings will be passed to
        the :class:`logging.Formatter` constructor.
    level : :class:`str`
        Logging level for the handler.

    Returns
    ----------
    :class:`logging.Logger`
        Configured logger.

    Notes
    -----
    Duplicated version of pyrolite.util.log.Handle.
    """
    if isinstance(logger, str):
        logger = logging.getLogger(logger)
    elif isinstance(logger, logging.Logger):
        pass
    else:
        raise NotImplementedError
    logger.propagate = False
    if isinstance(formatter, str):
        formatter = logging.Formatter(formatter)

    active_handlers = [
        i
        for i in logger.handlers
        if isinstance(i, (handler_class))  # not a null handler
    ]
    if active_handlers:
        handler = active_handlers[0]  # use the existing stream handler
    else:
        handler = handler_class()
    handler.setFormatter(formatter)
    if handler not in active_handlers:
        logger.addHandler(handler)
    if level is not None:
        logger.setLevel(getattr(logging, level))
    return logger
