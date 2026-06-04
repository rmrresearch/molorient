from decimal import Decimal, getcontext
from molorient.classes.atom import Atom
import numpy as np


def arccos_series(z):
    """
    Taylor series expansion of arccos around 0 with the equation
    arccos(z) = pi/2 - arcsin(z).
    """

    if z == Decimal('1'):
        return Decimal('0')
    
    elif z == Decimal('-1'):
        return Decimal('3.141592653589793238462643383')  # pi to 28 decimal places

    getcontext().prec += 2

    #Pi with 30 significant figures
    pi_dec = Decimal('3.141592653589793238462643383280')

    i, lasts, s, num, coeff = Decimal('0'), Decimal('0'), Decimal(z), Decimal(z), Decimal('1')
    tol = Decimal('1e-28')
    while True:
        lasts = s
        i += 1
        num *= z * z
        coeff *= ((2*i - 1)**2 / ((2*i) * (2*i + 1)))
        s += num * coeff
        if abs(s - lasts) < tol:
            break

    result = pi_dec / 2 - s
    getcontext().prec -= 2
    return +result


def cos_series(y):
    """
    Taylor series expansion of cos(y) around 0.
    """

    getcontext().prec += 2
    i, lasts, s, fact, num, sign = Decimal('0'), Decimal('0'), Decimal('1'), Decimal('1'), Decimal('1'), Decimal('1')
    tol = Decimal('1e-28')
    while True:
        lasts = s
        i += 2
        fact *= i * (i-1)
        num *= y * y
        sign *= -1
        s += num / fact * sign
        if abs(s - lasts) < tol:
            break

    getcontext().prec -= 2
    return +s


def inertia_tensor(atoms):
    """
    This function calculates the inertia tensor of the system using charge instead of mass. It returns
    the principal moments of inertia (eigenvalues) of the tensor in order of increasing magnitude.
    """
    I_xx = sum([atom.charge * (atom.y**2 + atom.z**2) for atom in atoms])
    I_yy = sum([atom.charge * (atom.x**2 + atom.z**2) for atom in atoms])
    I_zz = sum([atom.charge * (atom.x**2 + atom.y**2) for atom in atoms])
    I_xy = -sum([atom.charge * atom.x * atom.y for atom in atoms])
    I_xz = -sum([atom.charge * atom.x * atom.z for atom in atoms])
    I_yz = -sum([atom.charge * atom.y * atom.z for atom in atoms])

    # The cubic characteristic polynomial of the inertia tensor is given by:
    a = -1 
    b = I_xx + I_yy + I_zz
    c = I_xy**2 + I_xz**2 + I_yz**2 - I_xx*I_yy - I_xx*I_zz - I_yy*I_zz
    d = I_xx*I_yy*I_zz - 2*I_xy*I_xz*I_yz - I_xx*I_yz**2 - I_yy*I_xz**2 - I_zz*I_xy**2

    # Converting to a depressed cubic form: t^3 + pt + q = 0
    p = (3*a*c - b**2) / (3*a**2)
    q = (2*b**3 - 9*a*b*c + 27*a**2*d) / (27*a**3)