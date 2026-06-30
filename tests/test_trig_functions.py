import numpy as np
from molorient.classes.atom import Atom
from molorient.utils.trig_helpers import arccos_series, arcsin_series, cos_series, pi_as_decimal, sin_series, arctan_series, arctan2
from decimal import Decimal, getcontext, ROUND_HALF_UP


def test_pi_as_decimal():
    orig_prec = getcontext().prec
    getcontext().prec = 28

    tol = Decimal('1e-6')
    assert pi_as_decimal() - Decimal(np.pi) < tol

    getcontext().prec = orig_prec


def test_cos():
    orig_prec = getcontext().prec
    getcontext().prec = 28


    tol = Decimal('1e-27')
    # Test cos_series function with known values
    assert abs(cos_series(Decimal('0')) - Decimal('1')) < tol 
    assert abs(cos_series((pi_as_decimal())/2 - Decimal('0'))) < tol
    assert abs(cos_series(pi_as_decimal()) - Decimal('-1')) < tol
    
    getcontext().prec = orig_prec



def test_arcsin():
    orig_prec = getcontext().prec
    getcontext().prec = 28

    tol = Decimal('1e-27')
    # Test arcsin_series function with known values
    assert abs(arcsin_series(Decimal('0')) - Decimal('0')) < tol
    assert abs(arcsin_series(Decimal('0.5')) - pi_as_decimal() / 6) < tol
    assert abs(arcsin_series(Decimal('1')) - pi_as_decimal() / 2) < tol
    assert abs(arcsin_series(Decimal('-1')) + pi_as_decimal() / 2) < tol

    getcontext().prec = orig_prec


def test_arccos():
    original_prec = getcontext().prec

    # Test arccos_series function with an even range of values from -1 to 1
    # comparing values to truncated NumPy values.
    arr = np.linspace(-1, 1, num=5)
    for num in arr:
        dec_num = Decimal(str(num))
        expected = Decimal(str(np.arccos(num))).quantize(Decimal("0.0000001"), rounding=ROUND_HALF_UP)
        getcontext().prec = 8
        tol = Decimal(10) ** -Decimal(getcontext().prec)
        actual = arccos_series(dec_num)
        assert abs(actual - expected) < tol
        
    getcontext().prec = original_prec


def test_sin():
    orig_prec = getcontext().prec
    getcontext().prec = 28
    tol = Decimal(10) ** -Decimal(getcontext().prec - 2)

    assert abs(sin_series(Decimal('0') - Decimal('0'))) < tol
    assert abs(sin_series(pi_as_decimal() / 2 ) - Decimal('1')) < tol
    assert abs(sin_series(pi_as_decimal()) - Decimal('0')) < tol

    getcontext().prec = orig_prec


def test_arctan():
    orig_prec = getcontext().prec
    getcontext().prec = 28

    tol = Decimal(10) ** -Decimal(getcontext().prec - 2)

    assert abs(arctan_series(Decimal('0')) - Decimal('0')) < tol
    assert abs(arctan_series(Decimal('1')) - (pi_as_decimal()) / 4) < tol
    assert abs(arctan_series(Decimal('-1')) - (-pi_as_decimal()) / 4) < tol

    getcontext().prec = orig_prec


def test_arctan2():
    orig_prec = getcontext().prec
    getcontext().prec = 28

    tol = Decimal(10) ** -Decimal(getcontext().prec - 2)
    
    assert abs(arctan2(1, 1) - pi_as_decimal() / 4) < tol 
    assert abs(arctan2(1, -1) - 3 * pi_as_decimal() / 4) < tol
    assert abs(arctan2(0, Decimal(1)) - 0) < tol

    getcontext().prec = orig_prec