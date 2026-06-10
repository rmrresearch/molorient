from decimal import Decimal, getcontext, ROUND_HALF_UP
from molorient.classes.atom import Atom
from molorient.utils.trig_helpers import arccos_series, cos_series, pi_as_decimal, arcsin_series
import numpy as np


def diagonalize_3_by_3(row_1, row_2, row_3):
    """
    Diagonalizes a 3x3 Hermitian matrix using Viete's Trigonometric Method
    for the cubic characteristic polynomial.
    """

    getcontext().prec += 2

    e, f, g = row_1
    _, h, i = row_2
    _, _, j = row_3

    a = -1
    b = e + h + j
    c = f**2 + g**2 + i**2 - e*h - e*j - h*j
    d = e*h*j - 2*f*g*i - e*i**2 - h*g**2 - j*f**2

    # Converting to a depressed cubic form: t^3 + pt + q = 0
    # of which x = t - b/3a
    p = (3*a*c - b**2) / (3*a**2)
    q = (2*b**3 - 9*a*b*c + 27*a**2*d) / (27*a**3)

    #Debugging
    print("p: ", p, "q: ", q)

    z = (3 * q) / (2 * p) * (-3 / p).sqrt()
    y = arccos_series(z)

    #Debugging
    print("z: ", z, "y: ", y)

    t_0 = cos_series(y / 3)
    t_1 = cos_series((y / 3) - (2 * pi_as_decimal()) / 3)
    t_2 = cos_series((y / 3) - (4 * pi_as_decimal()) / 3)

    sqrt_term = (-p / 3).sqrt()
    coeff_term = -b / 3

    x_0 = (2 * sqrt_term * t_0) - coeff_term
    x_1 = (2 * sqrt_term * t_1) - coeff_term
    x_2 = (2 * sqrt_term * t_2) - coeff_term

    getcontext().prec -=2
    
    return x_0, x_1, x_2