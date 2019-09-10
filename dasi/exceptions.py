class DASiException(Exception):
    pass


class AlignmentException(DASiException):
    pass


class AlignmentContainerException(DASiException):
    pass


class DasiDesignException(DASiException):
    pass


from dasi.utils.npdf import NumpyDataFrameException
