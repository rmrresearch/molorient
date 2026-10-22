from decimal import Decimal, getcontext
import random
from molorient.classes.atom import Atom
from molorient.classes.vector import Vector
from molorient.classes.square_matrix import SquareMatrix
from molorient.orientation.standardize_axes.nwchem_orientation import orient_mol
from molorient.helpers.trig_helpers import cos_series, sin_series, pi_as_decimal


def translate_atoms(atoms, trans_vec):
    """
    Translates atoms by a fixed (Decimal) vector. Returns NEW Atom objects;
    does not modify the input list, since callers reuse the same atoms
    for multiple checks (rotation, translation, ...).
    """
    getcontext().prec += 5

    translated_atoms = [
        Atom(atom.element, atom.x + trans_vec[0], atom.y + trans_vec[1], atom.z + trans_vec[2], atom.charge)
        for atom in atoms
    ]

    getcontext().prec -= 5
    return translated_atoms


def random_translation_vector():
    """Random Decimal vector, +/-10 Angstrom, built from integers (no binary-float noise)."""
    return [Decimal(random.randint(-1000000, 1000000)) / Decimal(100000) for _ in range(3)]


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


def sorted_geometry(atoms):
    """Atoms in a permutation-free order (equivalent atoms may swap places)."""
    return sorted(atoms, key=lambda a: (a.element, round(a.x, 5), round(a.y, 5), round(a.z, 5)))


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
        for a, b in zip(sorted_geometry(oriented_atoms), sorted_geometry(ref_atoms))
    )

    rotated_atoms = orient_mol(rotate_atoms_small(atoms, 137, 61, 253))
    assert all(
        abs(a.x - b.x) < 1e-6 and abs(a.y - b.y) < 1e-6 and abs(a.z - b.z) < 1e-6
        for a, b in zip(sorted_geometry(rotated_atoms), sorted_geometry(ref_atoms))
    )

    reversed_atoms = orient_mol(rotate_atoms_small(atoms[::-1], 137, 61, 253))
    assert all(
        a.element == b.element and abs(a.x - b.x) < 1e-6 and abs(a.y - b.y) < 1e-6 and abs(a.z - b.z) < 1e-6
        for a, b in zip(sorted_geometry(reversed_atoms), sorted_geometry(ref_atoms))
    )

    translated_atoms = orient_mol(translate_atoms(atoms, random_translation_vector()))
    assert all(
        abs(a.x - b.x) < 1e-6 and abs(a.y - b.y) < 1e-6 and abs(a.z - b.z) < 1e-6
        for a, b in zip(sorted_geometry(translated_atoms), sorted_geometry(ref_atoms))
    )


def test_asymmetric_molecule():
    atoms =[
        Atom("O", 0, 0, Decimal('-0.1168520'), 8),
        Atom("H", 0, Decimal('0.7546680'), Decimal('0.4674100'), 1),
        Atom("H", 0, Decimal('-0.7546680'), Decimal('0.4674100'), 1)
    ]

    ref_atoms = [
        Atom("O", 0, 0, Decimal('0.11685240'), 8),
        Atom("H", Decimal('-0.75466800'), 0, Decimal('-0.46740960'), 1),
        Atom("H", Decimal('0.75466800'), 0, Decimal('-0.46740960'), 1)
    ]

    oriented_atoms = orient_mol(atoms)
    assert all(
        abs(a.x - b.x) < 1e-6 and abs(a.y - b.y) < 1e-6 and abs(a.z - b.z) < 1e-6
        for a, b in zip(sorted_geometry(oriented_atoms), sorted_geometry(ref_atoms))
    )

    rotated_atoms = orient_mol(rotate_atoms_small(atoms, 137, 61, 253))
    assert all(
        abs(a.x - b.x) < 1e-6 and abs(a.y - b.y) < 1e-6 and abs(a.z - b.z) < 1e-6
        for a, b in zip(sorted_geometry(rotated_atoms), sorted_geometry(ref_atoms))
    )

    reversed_atoms = orient_mol(rotate_atoms_small(atoms[::-1], 137, 61, 253))
    assert all(
        a.element == b.element and abs(a.x - b.x) < 1e-6 and abs(a.y - b.y) < 1e-6 and abs(a.z - b.z) < 1e-6
        for a, b in zip(sorted_geometry(reversed_atoms), sorted_geometry(ref_atoms))
    )

    translated_atoms = orient_mol(translate_atoms(atoms, random_translation_vector()))
    assert all(
        abs(a.x - b.x) < 1e-6 and abs(a.y - b.y) < 1e-6 and abs(a.z - b.z) < 1e-6
        for a, b in zip(sorted_geometry(translated_atoms), sorted_geometry(ref_atoms))
    )


def test_symmetric_molecule():
    atoms = [
        Atom("N", 0, Decimal(-0.1141360), 0, 7),
        Atom("H", 0, Decimal(0.2663170), Decimal(0.9363560), 1),
        Atom("H", Decimal(0.8109080), Decimal(0.2663170), Decimal(-0.4681780), 1),
        Atom("H", Decimal(-0.8109080), Decimal(0.2663170), Decimal(-0.4681780), 1)
    ]

    ref_atoms = [
        Atom("N", 0, 0, Decimal('0.11413590'), 7),
        Atom("H", Decimal('0.66210368'), Decimal('0.66210368'), Decimal('-0.26631710'), 1),
        Atom("H", Decimal('0.24234671'), Decimal('-0.90445038'), Decimal('-0.26631710'), 1),
        Atom("H", Decimal('-0.90445038'), Decimal('0.24234671'), Decimal('-0.26631710'), 1)
    ]

    oriented_atoms = orient_mol(atoms)
    assert all(
        abs(a.x - b.x) < 1e-6 and abs(a.y - b.y) < 1e-6 and abs(a.z - b.z) < 1e-6
        for a, b in zip(sorted_geometry(oriented_atoms), sorted_geometry(ref_atoms))
    )

    rotated_atoms = orient_mol(rotate_atoms_small(atoms, 137, 61, 253))
    assert all(
        abs(a.x - b.x) < 1e-6 and abs(a.y - b.y) < 1e-6 and abs(a.z - b.z) < 1e-6
        for a, b in zip(sorted_geometry(rotated_atoms), sorted_geometry(ref_atoms))
    )

    reversed_atoms = orient_mol(rotate_atoms_small(atoms[::-1], 137, 61, 253))
    assert all(
        a.element == b.element and abs(a.x - b.x) < 1e-6 and abs(a.y - b.y) < 1e-6 and abs(a.z - b.z) < 1e-6
        for a, b in zip(sorted_geometry(reversed_atoms), sorted_geometry(ref_atoms))
    )

    translated_atoms = orient_mol(translate_atoms(atoms, random_translation_vector()))
    assert all(
        abs(a.x - b.x) < 1e-6 and abs(a.y - b.y) < 1e-6 and abs(a.z - b.z) < 1e-6
        for a, b in zip(sorted_geometry(translated_atoms), sorted_geometry(ref_atoms))
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
        Atom("H", Decimal('0.62695100'), Decimal('-0.62695100'), Decimal('0.62695100'), 1),
        Atom("H", Decimal('-0.62695100'), Decimal('0.62695100'), Decimal('0.62695100'), 1),
        Atom("H", Decimal('-0.62695100'), Decimal('-0.62695100'), Decimal('-0.62695100'), 1),
        Atom("H", Decimal('0.62695100'), Decimal('0.62695100'), Decimal('-0.62695100'), 1)
    ]

    oriented_atoms = orient_mol(atoms)
    assert all(
        abs(a.x - b.x) < 1e-6 and abs(a.y - b.y) < 1e-6 and abs(a.z - b.z) < 1e-6
        for a, b in zip(sorted_geometry(oriented_atoms), sorted_geometry(ref_atoms))
    )

    rotated_atoms = orient_mol(rotate_atoms_small(atoms, 137, 61, 253))
    assert all(
        abs(a.x - b.x) < 1e-6 and abs(a.y - b.y) < 1e-6 and abs(a.z - b.z) < 1e-6
        for a, b in zip(sorted_geometry(rotated_atoms), sorted_geometry(ref_atoms))
    )

    reversed_atoms = orient_mol(rotate_atoms_small(atoms[::-1], 137, 61, 253))
    assert all(
        a.element == b.element and abs(a.x - b.x) < 1e-6 and abs(a.y - b.y) < 1e-6 and abs(a.z - b.z) < 1e-6
        for a, b in zip(sorted_geometry(reversed_atoms), sorted_geometry(ref_atoms))
    )

    translated_atoms = orient_mol(translate_atoms(atoms, random_translation_vector()))
    assert all(
        abs(a.x - b.x) < 1e-6 and abs(a.y - b.y) < 1e-6 and abs(a.z - b.z) < 1e-6
        for a, b in zip(sorted_geometry(translated_atoms), sorted_geometry(ref_atoms))
    )


def check_invariance(atoms, ref_atoms):
    """Same geometry after rotation, atom reversal and translation."""
    variants = [
        atoms,
        rotate_atoms_small(atoms, 137, 61, 253),
        rotate_atoms_small(atoms[::-1], 137, 61, 253),
        translate_atoms(atoms, random_translation_vector()),
    ]
    for variant in variants:
        assert all(
            a.element == b.element and abs(a.x - b.x) < 1e-6 and abs(a.y - b.y) < 1e-6 and abs(a.z - b.z) < 1e-6
            for a, b in zip(sorted_geometry(orient_mol(variant)), sorted_geometry(ref_atoms))
        )


def test_icosahedral_molecule():
    phi = (1 + Decimal(5).sqrt()) / 2
    s = Decimal('2.8') / (1 + phi**2).sqrt()
    atoms = [Atom("W", 0, 0, 0, 74)]
    for a in (1, -1):
        for b in (1, -1):
            atoms += [Atom("Au", 0, a * s, b * s * phi, 79),
                      Atom("Au", a * s, b * s * phi, 0, 79),
                      Atom("Au", b * s * phi, 0, a * s, 79)]

    ref_atoms = [Atom("W", 0, 0, 0, 74)] + [
        Atom("Au", Decimal(x), Decimal(y), Decimal(z), 79) for x, y, z in (
            ('-2.38182226', '1.47204711', '0'), ('-2.38182226', '-1.47204711', '0'),
            ('-1.47204711', '0', '-2.38182226'), ('0', '-2.38182226', '1.47204711'),
            ('1.47204711', '0', '-2.38182226'), ('0', '2.38182226', '1.47204711'),
            ('0', '2.38182226', '-1.47204711'), ('-1.47204711', '0', '2.38182226'),
            ('0', '-2.38182226', '-1.47204711'), ('2.38182226', '-1.47204711', '0'),
            ('2.38182226', '1.47204711', '0'), ('1.47204711', '0', '2.38182226'))]

    check_invariance(atoms, ref_atoms)


def test_heavy_octahedral_molecule():
    vertices = ((2, 0, 0), (-2, 0, 0), (0, 2, 0), (0, -2, 0), (0, 0, 2), (0, 0, -2))
    atoms = [Atom("U", 0, 0, 0, 92)] + [Atom("F", *(Decimal(c) for c in v), 9) for v in vertices]
    ref_atoms = [Atom("U", 0, 0, 0, 92)] + [Atom("F", *(Decimal(c) for c in v), 9) for v in vertices]

    check_invariance(atoms, ref_atoms)


def test_low_precision_does_not_raise():
    """canonical_key's fixed 1e-6 quantize needs headroom above the caller's precision."""
    orig_prec = getcontext().prec
    getcontext().prec = 6

    atoms = [
        Atom("O", 0, 0, Decimal('-0.1168520'), 8),
        Atom("H", 0, Decimal('0.7546680'), Decimal('0.4674100'), 1),
        Atom("H", 0, Decimal('-0.7546680'), Decimal('0.4674100'), 1)
    ]
    try:
        assert len(orient_mol(atoms)) == 3
    finally:
        getcontext().prec = orig_prec
