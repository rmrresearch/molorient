from decimal import Decimal, getcontext, ROUND_HALF_UP
from molorient.classes.atom import Atom
from molorient.utils.trig_helpers import arccos_series, cos_series, pi_as_decimal, arcsin_series
import numpy as np


def inertia_tensor(atoms):
    """
    This function calculates the inertia tensor of the system using charge instead of mass. It returns
    the principal moments of inertia (eigenvalues) of the tensor in order of increasing magnitude.
    """
    getcontext().prec += 2

    I_xx = sum([atom.charge * (atom.y**2 + atom.z**2) for atom in atoms])
    I_yy = sum([atom.charge * (atom.x**2 + atom.z**2) for atom in atoms])
    I_zz = sum([atom.charge * (atom.x**2 + atom.y**2) for atom in atoms])
    I_xy = -sum([atom.charge * atom.x * atom.y for atom in atoms])
    I_xz = -sum([atom.charge * atom.x * atom.z for atom in atoms])
    I_yz = -sum([atom.charge * atom.y * atom.z for atom in atoms])

    # row_1 = [I_xx, I_xy, I_xz]
    # row_2 = [I_xy, I_yy, I_yz]
    # row_3 = [I_xz, I_yz, I_zz]

    # x_0, x_1, x_2 = diagonalize_3_by_3(row_1, row_2, row_3)

    # The cubic characteristic polynomial of the inertia tensor is given by:
    a = -1 
    b = I_xx + I_yy + I_zz
    c = I_xy**2 + I_xz**2 + I_yz**2 - I_xx*I_yy - I_xx*I_zz - I_yy*I_zz
    d = I_xx*I_yy*I_zz - 2*I_xy*I_xz*I_yz - I_xx*I_yz**2 - I_yy*I_xz**2 - I_zz*I_xy**2

    # Converting to a depressed cubic form: t^3 + pt + q = 0
    # of which x = t - b/3a
    p = (3*a*c - b**2) / (3*a**2)
    q = (2*b**3 - 9*a*b*c + 27*a**2*d) / (27*a**3)
    print("q: ", q, "p: ", p)
    z = (3 * q) / (2 * p) * (-3 / p).sqrt()
    y = arccos_series(Decimal("-0.999999999999999999999"))
    print("y: ", y)

    t_0 = cos_series(y / 3)
    t_1 = cos_series((y / 3) - (2 * pi_as_decimal()) / 3)
    t_2 = cos_series((y / 3) - (4 * pi_as_decimal()) / 3)

    sqrt_term = (-p / 3).sqrt()
    coeff_term = -b / 3

    x_0 = (2 * sqrt_term * t_0) - coeff_term
    x_1 = (2 * sqrt_term * t_1) - coeff_term
    x_2 = (2 * sqrt_term * t_2) - coeff_term

    moments = []

    getcontext().prec -= 2
    tol = Decimal(10) ** -(Decimal(getcontext().prec - 1))

    moments = sorted([x_0.quantize(tol, rounding = ROUND_HALF_UP), x_1.quantize(tol, rounding = ROUND_HALF_UP), x_2.quantize(tol, rounding = ROUND_HALF_UP)])

    moment_a, moment_b, moment_c = moments

    print(moment_a, moment_b, moment_c)

    return moment_a, moment_b, moment_c