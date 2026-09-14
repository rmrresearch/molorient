from decimal import Decimal, getcontext
from molorient.classes.atom import Atom
from molorient.classes.vector import Vector
from molorient.classes.square_matrix import SquareMatrix
from molorient.utils.nwchem_orientation import orient_mol
from molorient.utils.trig_helpers import cos_series, sin_series, pi_as_decimal


def rotate_atoms_small(atoms, alpha_deg, beta_deg, gamma_deg):
    """
    Rotates atoms by a small, fixed Tait-Bryan rotation (angles in degrees).
    Uses the same rotation-matrix convention as test_random_rotations.py.
    """
    getcontext().prec += 10

    alpha = Decimal(alpha_deg) * pi_as_decimal() / 180
    beta = Decimal(beta_deg) * pi_as_decimal() / 180
    gamma = Decimal(gamma_deg) * pi_as_decimal() / 180

    rot_mat = SquareMatrix(3)
    rot_mat.assign(0, 0, cos_series(beta) * cos_series(gamma))
    rot_mat.assign(0, 1, -cos_series(beta) * sin_series(gamma))
    rot_mat.assign(0, 2, sin_series(beta))
    rot_mat.assign(1, 0, (cos_series(alpha) * sin_series(gamma)) + (sin_series(alpha) * sin_series(beta) * cos_series(gamma)))
    rot_mat.assign(1, 1, (cos_series(alpha) * cos_series(gamma)) - (sin_series(alpha) * sin_series(beta) * sin_series(gamma)))
    rot_mat.assign(1, 2, -sin_series(alpha) * cos_series(beta))
    rot_mat.assign(2, 0, (sin_series(alpha) * sin_series(gamma)) - (cos_series(alpha) * sin_series(beta) * cos_series(gamma)))
    rot_mat.assign(2, 1, (sin_series(alpha) * cos_series(gamma)) + (cos_series(alpha) * sin_series(beta) * sin_series(gamma)))
    rot_mat.assign(2, 2, cos_series(alpha) * cos_series(beta))

    rotated_atoms = []
    for atom in atoms:
        pos_vec = Vector(3)
        pos_vec.assign(0, atom.x)
        pos_vec.assign(1, atom.y)
        pos_vec.assign(2, atom.z)

        new_pos = rot_mat.transpose().multiply(pos_vec)

        rotated_atoms.append(Atom(
            atom.element, new_pos.elements[0], new_pos.elements[1], new_pos.elements[2], atom.charge
        ))

    getcontext().prec -= 10
    return rotated_atoms


def test_linear_molecule():
    atoms = [
        Atom("C", 0, 0, 0, 6),
        Atom("O", 0, Decimal('1.1604230'), 0, 8),
        Atom("O", 0, Decimal('-1.1604230'), 0, 8)
    ]

    ref_atoms = [
        Atom("C", 0, 0, 0, 6),
        Atom("O", 0, 0, Decimal('1.1604230'), 8),
        Atom("O", 0, 0, Decimal('-1.1604230'), 8)
    ]

    oriented_atoms = orient_mol(atoms)
    assert all(
        abs(a.x - b.x) < 1e-6 and abs(a.y - b.y) < 1e-6 and abs(a.z - b.z) < 1e-6
        for a, b in zip(oriented_atoms, ref_atoms)
    )

    rotated_atoms = orient_mol(rotate_atoms_small(atoms, 2, 2, 2))
    assert all(
        abs(a.x - b.x) < 1e-6 and abs(a.y - b.y) < 1e-6 and abs(a.z - b.z) < 1e-6
        for a, b in zip(rotated_atoms, ref_atoms)
    )


def test_asymmetric_molecule():
    atoms =[
        Atom("O", 0, 0, Decimal('-0.1168520'), 8),
        Atom("H", 0, Decimal('0.7546680'), Decimal('0.4674100'), 1),
        Atom("H", 0, Decimal('-0.7546680'), Decimal('0.4674100'), 1)
    ]

    ref_atoms = [
        Atom("O", 0, 0, Decimal('-0.1168520'), 8),
        Atom("H", Decimal(-0.7546680), 0, Decimal('0.4674100'), 1),
        Atom("H", Decimal(0.7546680), 0, Decimal('0.4674100'), 1)
    ]

    oriented_atoms = orient_mol(atoms)
    assert all(
        abs(a.x - b.x) < 1e-6 and abs(a.y - b.y) < 1e-6 and abs(a.z - b.z) < 1e-6
        for a, b in zip(oriented_atoms, ref_atoms)
    )

    rotated_atoms = orient_mol(rotate_atoms_small(atoms, 2, 2, 2))
    assert all(
        abs(a.x - b.x) < 1e-6 and abs(a.y - b.y) < 1e-6 and abs(a.z - b.z) < 1e-6
        for a, b in zip(rotated_atoms, ref_atoms)
    )


def test_symmetric_molecule():
    atoms = [
        Atom("N", 0, Decimal(-0.1141360), 0, 7),
        Atom("H", 0, Decimal(0.2663170), Decimal(0.9363560), 1),
        Atom("H", Decimal(0.8109080), Decimal(0.2663170), Decimal(-0.4681780), 1),
        Atom("H", Decimal(-0.8109080), Decimal(0.2663170), Decimal(-0.4681780), 1)
    ]

    ref_atoms = [
        Atom("N", 0, 0, Decimal('-0.11413590'), 7),
        Atom("H", Decimal('0.66210364'), Decimal('0.66210364'), Decimal('0.26631710'), 1),
        Atom("H", Decimal('-0.90445040'), Decimal('0.24234675'), Decimal('0.26631710'), 1),
        Atom("H", Decimal('0.24234675'), Decimal('-0.90445040'), Decimal('0.26631710'), 1)
    ]

    oriented_atoms = orient_mol(atoms)
    assert all(
        abs(a.x - b.x) < 1e-6 and abs(a.y - b.y) < 1e-6 and abs(a.z - b.z) < 1e-6
        for a, b in zip(oriented_atoms, ref_atoms)
    )

    rotated_atoms = orient_mol(rotate_atoms_small(atoms, 2, 2, 2))
    assert all(
        abs(a.x - b.x) < 1e-6 and abs(a.y - b.y) < 1e-6 and abs(a.z - b.z) < 1e-6
        for a, b in zip(rotated_atoms, ref_atoms)
    )


def test_spherical_molecule():
    atoms = [
        Atom("C", 0, 0, 0, 6),
        Atom("H", Decimal('0.6269510'), Decimal('0.6269510'), Decimal('-0.6269510'), 1),
        Atom("H", Decimal('-0.6269510'), Decimal('-0.6269510'), Decimal('-0.6269510'), 1),
        Atom("H", Decimal('-0.6269510'), Decimal('0.6269510'), Decimal('0.6269510'), 1),
        Atom("H", Decimal('0.6269510'), Decimal('-0.6269510'), Decimal('0.6269510'), 1)
    ]

    ref_atoms = [
        Atom("C", 0, 0, 0, 6),
        Atom("H", Decimal('0.62695100'), Decimal('-0.62695100'), Decimal('-0.62695100'), 1),
        Atom("H", Decimal('-0.62695100'), Decimal('0.62695100'), Decimal('-0.62695100'), 1),
        Atom("H", Decimal('0.62695100'), Decimal('0.62695100'), Decimal('0.62695100'), 1),
        Atom("H", Decimal('-0.62695100'), Decimal('-0.62695100'), Decimal('0.62695100'), 1)
    ]

    oriented_atoms = orient_mol(atoms)
    assert all(
        abs(a.x - b.x) < 1e-6 and abs(a.y - b.y) < 1e-6 and abs(a.z - b.z) < 1e-6
        for a, b in zip(oriented_atoms, ref_atoms)
    )

    rotated_atoms = orient_mol(rotate_atoms_small(atoms, 2, 2, 2))
    assert all(
        abs(a.x - b.x) < 1e-6 and abs(a.y - b.y) < 1e-6 and abs(a.z - b.z) < 1e-6
        for a, b in zip(rotated_atoms, ref_atoms)
    )