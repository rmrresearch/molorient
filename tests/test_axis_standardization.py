from molorient.classes.atom import Atom
from molorient.utils.axis_standardization import arccos_series, cos_series, inertia_tensor
from decimal import Decimal


def test_arccos():
    tol = Decimal('1e-27')
    # Test arccos_series function with known values
    assert abs(arccos_series(Decimal('0')) - Decimal('1.570796326794896619231321692')) < tol  # pi/2 to 28 decimal places
    assert abs(arccos_series(Decimal('0.5')) - Decimal('1.047197551196597746154214461')) < tol  # arccos(0.5) = pi/3 to 28 significant figures
    assert abs(arccos_series(Decimal('1')) - Decimal('0')) < tol  # arccos(1) = 0


def test_cos():
    tol = Decimal('1e-27')
    # Test cos_series function with known values
    assert abs(cos_series(Decimal('0')) - Decimal('1')) < tol  # cos(0) = 1
    assert abs(cos_series(Decimal('1.570796326794896619231321692')) - Decimal('0')) < tol  # cos(pi/2) = 0 to 28 decimal places
    assert abs(cos_series(Decimal('3.141592653589793238462643383')) - Decimal('-1')) < tol  # cos(pi) = -1 to 28 decimal places