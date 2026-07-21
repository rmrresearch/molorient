from molorient.utils.axis_standardization import inertia_tensor, standardize_axes, cn_axes_finder
from molorient.classes.atom import Atom
import numpy as np
from decimal import Decimal, getcontext


def test_inertia_tensor():
    atoms =[
        Atom("H", 0.0, 0.0, 0.0, 1.0),
        Atom("H", 1.0, 0.0, 0.0, 1.0)
    ]

    moments, eigvecs = inertia_tensor(atoms) 

    vecs = np.array([
        [1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.0, 1.0 ,0.0]
    ])

    tol = Decimal('10')**(Decimal('-6'))
    assert (moments[0] - 0) < tol
    assert (moments[1] - 1) < tol
    assert (moments[2] - 1) < tol
    for i in range(3):
        assert eigvecs[0].elements[i] - Decimal(vecs[0][i]) < tol
        assert eigvecs[1].elements[i] - Decimal(vecs[1][i]) < tol
        assert eigvecs[2].elements[i] - Decimal(vecs[2][i]) < tol


def test_single_atom():
    atom = [
        Atom("H", 0.0, 0.0, 0.0, 1.0)
    ]

    moments, eigvecs = inertia_tensor(atom)

    assert all(m == 0 for m in moments)
    assert eigvecs[0].elements[0] == Decimal('1')
    assert eigvecs[1].elements[1] == Decimal('1')
    assert eigvecs[2].elements[2] == Decimal('1')

    std_atoms = standardize_axes(moments, eigvecs, atom)

    assert std_atoms[0].x == 0
    assert std_atoms[0].y == 0
    assert std_atoms[0].z == 0


def test_asymmetric_top():
    orig_prec = getcontext().prec
    getcontext().prec = 28
    atoms = [
        Atom("H", 0.0, 0.0, 1.0, 1.0),
        Atom("H", 0.0, 1.0, 0.0, 1.0),
        Atom("H", 1.0, 0.0, 1.0, 1.0),
    ]

    moments, eigvecs = inertia_tensor(atoms)
    std_atoms = standardize_axes(moments, eigvecs, atoms)
        
    expected_coords = [
        (Decimal('-1.3763819204711735382072095819'), 0, Decimal('-0.3249196962329063261558714122')),
        (Decimal('-0.8506508083520399321815404971'), 0, Decimal('0.5257311121191336060256690848')),
        (0, -1, Decimal(0))
    ]

    for atom in std_atoms:
        print(atom.element, atom.x, atom.y, atom.z)

    for atom, (x_exp, y_exp, z_exp) in zip(std_atoms, expected_coords):
        assert atom.x == x_exp
        assert atom.y == y_exp
        assert atom.z == z_exp
    
    getcontext().prec = orig_prec


def test_linear():
    atoms = [
        Atom("H", 0.0, 0.0, 0.0, 1.0),
        Atom("H", 0.0, 1.0, 0.0, 1.0),
    ]

    moments, eigvecs = inertia_tensor(atoms)
    std_atoms = standardize_axes(moments, eigvecs, atoms)

    expected_coords = [
        (0, 0, 0),
        (0, 0, -1)
    ]

    for atom, (x_exp, y_exp, z_exp) in zip(std_atoms, expected_coords):
        assert atom.x == x_exp
        assert atom.y == y_exp
        assert atom.z == z_exp


def test_symmetric():
    orig_prec = getcontext().prec
    getcontext().prec = 28
    atoms = [
        Atom("He", 0, 1, 0, 2),
        Atom("H", 0, -1, 1, 1),
        Atom("H", -Decimal(3).sqrt() / 2, -1, Decimal('-0.5'), 1),
        Atom("H", Decimal(3).sqrt() / 2, -1, Decimal('-0.5'), 1),
    ]

    moments, eigvecs = inertia_tensor(atoms)
    std_atoms = standardize_axes(moments, eigvecs, atoms)

    for atom in std_atoms:
        print(atom.x, atom.y, atom.z)

    expected_coords = [
        (0, 0, -1),
        (0, 1, 1),
        (-Decimal(3).sqrt() / 2, -Decimal('0.5'), 1),
        (Decimal(3).sqrt() / 2, -Decimal('0.5'), 1)
    ]

    for atom, (x_exp, y_exp, z_exp) in zip(std_atoms, expected_coords):
        assert atom.x == x_exp
        assert atom.y == y_exp
        assert atom.z == z_exp
    
    getcontext().prec = orig_prec


def test_spherical():
    orig_prec = getcontext().prec
    getcontext().prec = 28
    atoms = [
        Atom("C", 0, 0, 0, 6),
        Atom("H", Decimal('0.62911'), Decimal('0.62911'), Decimal('0.62911'), 1),
        Atom("H", Decimal('-0.62911'), Decimal('-0.62911'), Decimal('0.62911'), 1),
        Atom("H", Decimal('-0.62911'), Decimal('0.62911'), Decimal('-0.62911'), 1),
        Atom("H", Decimal('0.62911'), Decimal('-0.62911'), -Decimal('0.62911'), 1)
    ]

    moments, eigvecs = inertia_tensor(atoms)
    std_atoms = standardize_axes(moments, eigvecs, atoms) 
    for a in std_atoms:
        print(a.element, a.x, a.y, a.z)

    for i in range(5):
        assert atoms[i].x == std_atoms[i].x
        assert atoms[i].y == std_atoms[i].y
        assert atoms[i].z == std_atoms[i].z

    getcontext().prec = orig_prec


def test_cn_axes_finder():
    atoms = [
        Atom("N", 0, 0, 0, 7),
        Atom("H", Decimal('0.5939'), Decimal('0.5939'), Decimal('0.5939'), 7),
        Atom("H", -Decimal('0.5939'), -Decimal('0.5939'), Decimal('0.5939'), 1),
        Atom("H", -Decimal('0.5939'), Decimal('0.5939'), -Decimal('0.5939'), 1),
        Atom("H", Decimal('0.5939'), -Decimal('0.5939'), -Decimal('0.5939'), 1)
    ]

    group, c3 = cn_axes_finder(atoms)
    assert group == 'Td'
    assert c3[0].elements[2] == 1
    assert c3[1].elements[1] == 1
    assert c3[2].elements[0] == 1
