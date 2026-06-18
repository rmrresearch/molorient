from decimal import Decimal, getcontext
from molorient.classes.atom import Atom
from molorient.classes.square_matrix import SquareMatrix
from molorient.utils.diagonalization import eigval_solver, eigvec_solver
import numpy as np


def inertia_tensor(atoms):
    """
    This function calculates the inertia tensor of the system using charge instead of mass. It returns
    the principal moments of inertia (eigenvalues) of the tensor in order of increasing magnitude.
    """

    tensor = SquareMatrix(3)
    I_xx = sum([atom.charge * (atom.y**2 + atom.z**2) for atom in atoms])
    I_yy = sum([atom.charge * (atom.x**2 + atom.z**2) for atom in atoms])
    I_zz = sum([atom.charge * (atom.x**2 + atom.y**2) for atom in atoms])
    I_xy = -sum([atom.charge * atom.x * atom.y for atom in atoms])
    I_xz = -sum([atom.charge * atom.x * atom.z for atom in atoms])
    I_yz = -sum([atom.charge * atom.y * atom.z for atom in atoms])

    tensor.assign(0, 0, I_xx)
    tensor.assign(0, 1, I_xy)
    tensor.assign(0, 2, I_xz)
    tensor.assign(1, 0, I_xy)
    tensor.assign(1, 1, I_yy)
    tensor.assign(1, 2, I_yz)
    tensor.assign(2, 0, I_xz)
    tensor.assign(2, 1, I_yz)
    tensor.assign(2, 2, I_zz)

    moment_a, moment_b, moment_c = sorted(eigval_solver(tensor))
    v_0, v_1, v_2 = eigvec_solver(moment_a, moment_b, moment_c, tensor)
    eigvals = [moment_a, moment_b, moment_c]
    eigvecs = [v_0, v_1, v_2]
    print("eigvals: ", eigvals)
    print("v_0: ", eigvecs[0].elements)
    print("v_1: ", eigvecs[1].elements)
    print("v_2: ", eigvecs[2].elements)
    return eigvals, eigvecs