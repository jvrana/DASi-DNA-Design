import logging


logging.basicConfig(format="%(levelname)s: %(message)s", level=logging.WARNING)
logger = logging.getLogger("DASi")
logger.setLevel(logging.DEBUG)
