"""Logging configuration for Brimer-PLAST.

Sets up a logger that prints structured progress to stderr, controlled
by a ``--verbose`` / ``-v`` flag (repeatable).

Level scheme:

* ``-v`` (INFO):  pipeline-stage announcements (reading genome, designing, filtering)
* ``-vv`` (DEBUG): per-pair data (fragment lists, template coordinates)
"""

from __future__ import annotations

import logging
import sys

_LOG_CONFIGURED = False


def configure_logging(verbosity: int = 0) -> None:
    """Configure the ``brimer_plast`` logger.

    Args:
        verbosity: 0 = warning+error, 1 = info, 2 = debug.
    """
    global _LOG_CONFIGURED
    if _LOG_CONFIGURED:
        return
    _LOG_CONFIGURED = True

    level = logging.WARNING
    if verbosity >= 2:
        level = logging.DEBUG
    elif verbosity >= 1:
        level = logging.INFO

    logger = logging.getLogger("brimer_plast")
    logger.setLevel(level)

    handler = logging.StreamHandler(sys.stderr)
    handler.setLevel(level)
    handler.setFormatter(logging.Formatter("%(levelname)s: %(message)s"))

    # Remove any pre-existing handlers (e.g. from pytest)
    logger.handlers.clear()
    logger.addHandler(handler)


def get_logger() -> logging.Logger:
    """Return the ``brimer_plast`` logger, configuring defaults if needed."""
    configure_logging()
    return logging.getLogger("brimer_plast")
