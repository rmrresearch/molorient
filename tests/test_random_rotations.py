from molorient.utils.cli import parse_xyz, set_precision
from molorient.utils.orient_system import orient_system
from molorient.classes.square_matrix import SquareMatrix
from molorient.classes.vector import Vector
from molorient.classes.atom import Atom
from molorient.utils.trig_helpers import cos_series, sin_series, pi_as_decimal
from unittest.mock import patch
from decimal import Decimal, getcontext
import random


def test_rotate_asymmetric():
    """
    Rotates acetic_acid_ref.xyz about randomly generated Tait-Bryan angles.
    """

    atoms, folder, base, ext = parse_xyz("tests/randomized_tests/rotations/acetic_acid_ref.xyz")
    oriented_ref_atoms, folder, base, ext = parse_xyz("tests/randomized_tests/rotations/acetic_acid_oriented_6.xyz")
    with patch('builtins.input', return_value = '6'):
        result = set_precision()
        assert result == 6
    for atom in atoms:
        atom.x = +atom.x
        atom.y = +atom.y
        atom.z = +atom.z
    getcontext().prec += 10

    rot_mat = SquareMatrix(3)
    alpha = Decimal(random.random()) * (Decimal(2) * pi_as_decimal())
    beta = Decimal(random.random()) * pi_as_decimal()
    gamma = Decimal(random.random()) * (Decimal(2) * pi_as_decimal())
    rot_mat.assign(0, 0, cos_series(beta) * cos_series(gamma))
    rot_mat.assign(0, 1, -cos_series(beta) * sin_series(gamma))
    rot_mat.assign(0, 2, sin_series(beta))
    rot_mat.assign(1, 0, (cos_series(alpha) * sin_series(gamma)) + (sin_series(alpha) * sin_series(beta) * cos_series(gamma)))
    rot_mat.assign(1, 1, (cos_series(alpha) * cos_series(gamma)) - (sin_series(alpha) * sin_series(beta) * sin_series(gamma)))
    rot_mat.assign(1, 2, -sin_series(alpha) * cos_series(beta))
    rot_mat.assign(2, 0, (sin_series(alpha) * sin_series(gamma)) - (cos_series(alpha) * sin_series(beta) * cos_series(gamma)))
    rot_mat.assign(2, 1, (sin_series(alpha) * cos_series(gamma)) + (cos_series(alpha) * sin_series(beta) * sin_series(gamma)))
    rot_mat.assign(2, 2, cos_series(alpha) * cos_series(beta))

    #Rotation
    rotated_atoms = []

    for atom in atoms:
        pos_vec = Vector(3)
        pos_vec.assign(0, atom.x)
        pos_vec.assign(1, atom.y) 
        pos_vec.assign(2, atom.z)

        new_pos = (rot_mat.transpose()).multiply(pos_vec)

        rotated_atoms.append(Atom(atom.element, 
                                        new_pos.elements[0],
                                        new_pos.elements[1],
                                        new_pos.elements[2],
                                        atom.charge))
    
    getcontext().prec -= 10
    std_atoms = orient_system(rotated_atoms)

    for atom, ref in zip(std_atoms, oriented_ref_atoms):
        assert atom.element == ref.element
        assert atom.x == ref.x
        assert atom.y == ref.y
        assert atom.z == ref.z


def test_rotate_symmetric():
    """
    Rotates allene_ref.xyz about randomly generated Tait-Bryan angles
    """
    atoms, folder, base, ext = parse_xyz("tests/randomized_tests/rotations/allene_ref.xyz")
    oriented_ref_atoms, folder, base, ext = parse_xyz("tests/randomized_tests/rotations/allene_oriented_6.xyz")

    for atom in atoms:
        atom.x = +atom.x
        atom.y = +atom.y
        atom.z = +atom.z

    with patch('builtins.input', return_value = '6'):
        result = set_precision()
        assert result == 6

    getcontext().prec += 10

    rot_mat = SquareMatrix(3)
    alpha = Decimal(random.random()) * (Decimal(2) * pi_as_decimal())
    beta = Decimal(random.random()) * pi_as_decimal()
    gamma = Decimal(random.random()) * (Decimal(2) * pi_as_decimal())
    rot_mat.assign(0, 0, cos_series(beta) * cos_series(gamma))
    rot_mat.assign(0, 1, -cos_series(beta) * sin_series(gamma))
    rot_mat.assign(0, 2, sin_series(beta))
    rot_mat.assign(1, 0, (cos_series(alpha) * sin_series(gamma)) + (sin_series(alpha) * sin_series(beta) * cos_series(gamma)))
    rot_mat.assign(1, 1, (cos_series(alpha) * cos_series(gamma)) - (sin_series(alpha) * sin_series(beta) * sin_series(gamma)))
    rot_mat.assign(1, 2, -sin_series(alpha) * cos_series(beta))
    rot_mat.assign(2, 0, (sin_series(alpha) * sin_series(gamma)) - (cos_series(alpha) * sin_series(beta) * cos_series(gamma)))
    rot_mat.assign(2, 1, (sin_series(alpha) * cos_series(gamma)) + (cos_series(alpha) * sin_series(beta) * sin_series(gamma)))
    rot_mat.assign(2, 2, cos_series(alpha) * cos_series(beta))

    #Rotation
    rotated_atoms = []

    for atom in atoms:
        pos_vec = Vector(3)
        pos_vec.assign(0, atom.x)
        pos_vec.assign(1, atom.y) 
        pos_vec.assign(2, atom.z)

        new_pos = (rot_mat.transpose()).multiply(pos_vec)

        rotated_atoms.append(Atom(atom.element, 
                                        new_pos.elements[0],
                                        new_pos.elements[1],
                                        new_pos.elements[2],
                                        atom.charge))
    
    getcontext().prec -= 10
    std_atoms = orient_system(rotated_atoms)

    for atom, ref in zip(std_atoms, oriented_ref_atoms):
        assert atom.element == ref.element
        assert atom.x == ref.x
        assert atom.y == ref.y
        assert atom.z == ref.z


def test_rotate_spherical():
    """
    Rotates si_c4_h12.xyz about randomly generated Tait-Bryan angles.
    """
    atoms, folder, base, ext = parse_xyz("tests/randomized_tests/rotations/si_c4_h12.xyz")
    oriented_ref_atoms, folder, base, ext = parse_xyz("tests/randomized_tests/rotations/si_c4_h12_oriented_6.xyz")
    with patch('builtins.input', return_value = '6'):
        result = set_precision()
        assert result == 6
    
    for atom in atoms:
        atom.x = +atom.x
        atom.y = +atom.y
        atom.z = +atom.z
    
    getcontext().prec += 10

    rot_mat = SquareMatrix(3)
    alpha = Decimal(random.random()) * (Decimal(2) * pi_as_decimal())
    beta = Decimal(random.random()) * pi_as_decimal()
    gamma = Decimal(random.random()) * (Decimal(2) * pi_as_decimal())
    rot_mat.assign(0, 0, cos_series(beta) * cos_series(gamma))
    rot_mat.assign(0, 1, -cos_series(beta) * sin_series(gamma))
    rot_mat.assign(0, 2, sin_series(beta))
    rot_mat.assign(1, 0, (cos_series(alpha) * sin_series(gamma)) + (sin_series(alpha) * sin_series(beta) * cos_series(gamma)))
    rot_mat.assign(1, 1, (cos_series(alpha) * cos_series(gamma)) - (sin_series(alpha) * sin_series(beta) * sin_series(gamma)))
    rot_mat.assign(1, 2, -sin_series(alpha) * cos_series(beta))
    rot_mat.assign(2, 0, (sin_series(alpha) * sin_series(gamma)) - (cos_series(alpha) * sin_series(beta) * cos_series(gamma)))
    rot_mat.assign(2, 1, (sin_series(alpha) * cos_series(gamma)) + (cos_series(alpha) * sin_series(beta) * sin_series(gamma)))
    rot_mat.assign(2, 2, cos_series(alpha) * cos_series(beta))

    #Rotation
    rotated_atoms = []

    for atom in atoms:
        pos_vec = Vector(3)
        pos_vec.assign(0, atom.x)
        pos_vec.assign(1, atom.y) 
        pos_vec.assign(2, atom.z)

        new_pos = (rot_mat.transpose()).multiply(pos_vec)

        rotated_atoms.append(Atom(atom.element, 
                                        new_pos.elements[0],
                                        new_pos.elements[1],
                                        new_pos.elements[2],
                                        atom.charge))
    
    getcontext().prec -= 10
    std_atoms = orient_system(rotated_atoms)

    for atom, ref in zip(std_atoms, oriented_ref_atoms):
        assert atom.element == ref.element
        assert atom.x == ref.x
        assert atom.y == ref.y
        assert atom.z == ref.z