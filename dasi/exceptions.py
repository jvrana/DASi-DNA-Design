"""DASi exceptions and warnings"""

class DASiWarning(Warning):
    """A generic dasi warning"""


class DASiException(Exception):
    """A generic dasi exception"""


class AlignmentException(DASiException):
    """An alignment exception"""


class AlignmentContainerException(DASiException):
    """An alignment container exception"""


class DasiDesignException(DASiException):
    """A design exception"""
