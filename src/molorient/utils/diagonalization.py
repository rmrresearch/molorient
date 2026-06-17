from decimal import Decimal, getcontext, ROUND_HALF_UP
from molorient.utils.trig_helpers import arccos_series, cos_series, pi_as_decimal
from molorient.classes.vector import Vector
from molorient.classes.square_matrix import SquareMatrix
import numpy as np


def eigval_solver(squarematrix):
    """
    Finds eigenvalues of a 3x3 Hermitian Matrix using Viète's Trigonometric
    Method for the characteristic polynomial.
    """

    getcontext().prec += 2

    e, f, g = squarematrix.elements[0]
    _, h, i = squarematrix.elements[1]
    _, _, j = squarematrix.elements[2]

    a = -1
    b = e + h + j
    c = f**2 + g**2 + i**2 - e*h - e*j - h*j
    d = e*h*j + 2*f*g*i - e*i**2 - h*g**2 - j*f**2

    # Converting to a depressed cubic form: t^3 + pt + q = 0
    # of which x = t - b/3a
    p = (3*a*c - b**2) / (3*a**2)
    q = (2*b**3 - 9*a*b*c + 27*a**2*d) / (27*a**3)

    z = (3 * q) / (2 * p) * (-3 / p).sqrt()
    y = arccos_series(z)

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


def eigvec_solver(e_0, e_1, e_2, squarematrix):
    """
    Solves for eigenvectors of 3x3 Matrix.
    """