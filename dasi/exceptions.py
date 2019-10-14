"""DASi exceptions and warnings."""
from .log import logger
from dasi.utils.networkx.exceptions import NetworkxUtilsException


class DASiWarning(Warning):
    """A generic dasi warning."""


class DASiException(Exception):
    """A generic dasi exception."""


class AlignmentException(DASiException):
    """An alignment exception."""


class AlignmentContainerException(DASiException):
    """An alignment container exception."""


class DasiDesignException(DASiException):
    """A design exception."""


class DasiSequenceDesignException(DASiException):
    """Sequence design exception."""


class DasiNoPrimerPairsException(DASiException):
    """Sequence design exception."""
