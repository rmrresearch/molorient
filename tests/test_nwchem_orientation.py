from decimal import Decimal
from molorient.classes.atom import Atom
from molorient.utils.nwchem_orientation import orient_mol


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