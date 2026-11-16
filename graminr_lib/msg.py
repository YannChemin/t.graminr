"""
MODULE:     t.graminr (graminr_lib.msg)

AUTHOR(S):  Yann Chemin

PURPOSE:    User messaging that works both inside a GRASS session and out of
            one. Inside a session everything goes through grass.script, so
            messages honour --quiet/--verbose and MetaMonitor sees them;
            outside one (a pytest run of the model code, for instance) the
            same calls fall back to stderr and plain exceptions, which keeps
            the ported models testable without a mapset.

SPDX-License-Identifier: Unlicense
"""

import sys

try:
    import grass.script as gs
except ImportError:  # pragma: no cover - exercised only outside GRASS
    gs = None


class GraminrError(Exception):
    """Raised instead of calling G_fatal_error() outside a GRASS session."""


def fatal(message):
    """Abort loudly: gs.fatal() in a session, an exception outside one."""
    if gs is not None:
        gs.fatal(message)
    raise GraminrError(message)


def warning(message):
    if gs is not None:
        gs.warning(message)
    else:
        sys.stderr.write("WARNING: {}\n".format(message))


def message(text):
    if gs is not None:
        gs.message(text)
    else:
        sys.stderr.write("{}\n".format(text))


def verbose(text):
    if gs is not None:
        gs.verbose(text)


def percent(current, total, step=1):
    if gs is not None:
        gs.percent(current, total, step)


def translate(text):
    """gettext when GRASS installed it, the text itself otherwise."""
    try:
        return _(text)  # noqa: F821 - installed into builtins by grass.script
    except NameError:
        return text
