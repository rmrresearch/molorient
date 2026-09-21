from decimal import Decimal, getcontext


def prec_tol(offset=0):
    """
    10**-(current Decimal precision - offset): the tolerance / quantize target
    used throughout this codebase for "how close counts as equal", scaled to
    whatever Decimal precision is active when this is called. offset controls
    how many digits of margin are kept below full precision.
    """
    return Decimal(10) ** -(getcontext().prec - offset)
