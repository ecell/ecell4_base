from logging import getLogger
_log = getLogger(__name__)

try:
    import pint
except ImportError:
    HAS_PINT = False
    _log.warn("No module named 'pint' required by '{}'".format(__name__))
else:
    HAS_PINT = True
    from ._unit import *
