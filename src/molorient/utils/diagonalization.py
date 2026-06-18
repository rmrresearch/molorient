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

    getcontext().prec -= 2
    
    return x_0, x_1, x_2


def eigvec_solver(eig_0, eig_1, eig_2, squarematrix):
    """
    Solves for eigenvectors of 3x3 Hermitian matrix using the cross-product method to
    cut down on arithmetic operations.
    """

    getcontext().prec += 2

    a = squarematrix
    v_0 = Vector(3)
    v_1 = Vector(3)
    v_2 = Vector(3)
    a_0 = Vector(3)
    a_1 = Vector(3)
    a_2 = Vector(3)
    e_0 = Vector(3)
    e_1 = Vector(3)
    id_mat = SquareMatrix(3)
    cross_term_0 = Vector(3)
    cross_term_1 = Vector(3)

    e_0.assign(0, Decimal('1.0'))
    e_1.assign(1, Decimal('1.0'))
    for i in range(3):
        id_mat.elements[i][i] = Decimal('1.0')

    for i in range(3):
        a_0.elements[i] = a.elements[i][0]
        a_1.elements[i] = a.elements[i][1]
        a_2.elements[i] = a.elements[i][2]

        cross_term_0.elements[i] = a_0.elements[i] + (e_0.scale(-eig_0)).elements[i]
        cross_term_1.elements[i] = a_1.elements[i] + (e_1.scale(-eig_0)).elements[i]

    if all(x == 0 for x in cross_term_0.elements) or all(x == 0 for x in cross_term_1.elements):
            v_0.assign(0, Decimal('1'))
        
    mus = []

    for i in range(3):
        if cross_term_1.elements[i] != Decimal('0'):
            mu = cross_term_0.elements[i] / cross_term_1.elements[i]
            mus.append(mu)
    
    if all(m == mus[0] for m in mus):
        scale_term = Decimal('1.0') / (Decimal('1.0') + mus[0]**2).sqrt()
        w = Vector(3)
        w.assign(0, Decimal('1.0'))
        w.assign(1, -mus[0])
        v_0 = w.scale(scale_term)
    
    else:
        v_0 = cross_term_0.cross(cross_term_1)
    
    if eig_0 == eig_1:
        char_mat = SquareMatrix(3)
        char_mat_0 = Vector(3)

        for i in range(3):
            for j in range(3):
                char_mat.elements[i][j] = a.elements[i][j] + (id_mat.scale(-eig_0)).elements[i][j]
            char_mat_0.elements[i] = char_mat.elements[i][0]
        v_1 = v_0.cross(char_mat_0)
        
    else:
        for i in range(3):
            cross_term_0.elements[i] = a_0.elements[i] + (e_0.scale(-eig_1)).elements[i]
            cross_term_1.elements[i] = a_1.elements[i] + (e_1.scale(-eig_1).elements[i])

        v_1 = cross_term_0.cross(cross_term_1)

    v_2 = v_0.cross(v_1)

    vecs = [v_0, v_1, v_2]
    norm_vecs = []
    for v in vecs:
        norm = Decimal('1') / (v.elements[0]**2 + v.elements[1]**2 + v.elements[2]**2).sqrt()
        norm_v = Vector(3)
        for i in range(3):
            norm_v.elements[i] = (v.scale(norm)).elements[i]
        norm_vecs.append(norm_v)
    
    getcontext().prec -= 2

    return norm_vecs