from molorient.utils.precision import prec_tol
from decimal import Decimal, getcontext


def test_prec_tol():
    orig_prec = getcontext().prec
    getcontext().prec = 28

    assert prec_tol() == Decimal(10) ** -28
    assert prec_tol(2) == Decimal(10) ** -26

    getcontext().prec = orig_prec
