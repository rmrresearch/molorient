from molorient.orientation.standardize_axes.axis_standardization import inertia_tensor, standardize_axes, cn_axes_finder
from molorient.classes.atom import Atom
from decimal import Decimal, getcontext


def test_inertia_tensor():
    atoms =[
        Atom("H", 0.0, 0.0, 0.0, 1.0),
        Atom("H", 1.0, 0.0, 0.0, 1.0)
    ]

    moments, eigvecs = inertia_tensor(atoms) 

    tol = Decimal('10')**(Decimal('-6'))
    assert abs(moments[0] - 0) < tol
    assert abs(moments[1] - 1) < tol
    assert abs(moments[2] - 1) < tol

    #The molecule's axis is the (sign-arbitrary) eigenvector of the zero moment.
    assert abs(abs(eigvecs[0].elements[0]) - 1) < tol
    assert abs(eigvecs[0].elements[1]) < tol
    assert abs(eigvecs[0].elements[2]) < tol

    #The degenerate pair is an arbitrary orthonormal, right-handed completion.
    for i in range(3):
        assert abs(eigvecs[i].dot(eigvecs[i]) - 1) < tol
        for j in range(i):
            assert abs(eigvecs[i].dot(eigvecs[j])) < tol
    assert (eigvecs[0].cross(eigvecs[1])).dot(eigvecs[2]) > 0


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

    #The Td axes are the Cartesian axes, but which of the equivalent frames is used is
    #fixed by pick_canonical_frame (it may be the 90-degree-turned tetrahedron).
    assert std_atoms[0].element == "C" and (std_atoms[0].x, std_atoms[0].y, std_atoms[0].z) == (0, 0, 0)
    hydrogens = [a for a in std_atoms if a.element == "H"]
    assert len(hydrogens) == 4
    for h in hydrogens:
        for c in (h.x, h.y, h.z):
            assert abs(abs(c) - Decimal('0.62911')) < Decimal('1e-20')
    assert len({(h.x > 0, h.y > 0, h.z > 0) for h in hydrogens}) == 4
    assert len({h.x * h.y * h.z > 0 for h in hydrogens}) == 1   #one tetrahedron, not a cube corner mix

    getcontext().prec = orig_prec


def test_cn_axes_finder():
    atoms = [
        Atom("N", 0, 0, 0, 7),
        Atom("H", Decimal('0.5939'), Decimal('0.5939'), Decimal('0.5939'), 1),
        Atom("H", -Decimal('0.5939'), -Decimal('0.5939'), Decimal('0.5939'), 1),
        Atom("H", -Decimal('0.5939'), Decimal('0.5939'), -Decimal('0.5939'), 1),
        Atom("H", Decimal('0.5939'), -Decimal('0.5939'), -Decimal('0.5939'), 1)
    ]

    group, c3 = cn_axes_finder(atoms)
    assert group == 'Td'
    assert c3[0].elements[2] == 1
    assert c3[1].elements[1] == 1
    assert c3[2].elements[0] == 1


def test_classify_moments_uses_tolerance_tied_to_precision():
    from molorient.orientation.standardize_axes.axis_standardization import classify_moments
    D = Decimal
    assert classify_moments([D(0), D(0), D(0)], 6) == ('atom', None)
    assert classify_moments([D('1e-20'), D(50), D(50) + D('1e-6')], 6) == ('linear', 0)
    assert classify_moments([D(10), D(10) + D('1e-4') * 5, D(20)], 6) == ('symmetric', 2)
    assert classify_moments([D(10), D(20), D(20) + D('1e-6')], 6) == ('symmetric', 0)
    assert classify_moments([D(10), D(10) + D('1e-6'), D(10) + D('2e-6')], 6) == ('spherical', None)
    # just outside the tolerance (relative 1e-4 at p = 6) is asymmetric
    assert classify_moments([D(10), D(20), D(20) + D('0.01')], 6)[0] == 'asymmetric'
    # the same gap is a symmetric top at lower precision
    assert classify_moments([D(10), D(20), D(20) + D('0.01')], 4) == ('symmetric', 0)


def test_pick_canonical_frame_rejects_improper_frame_and_ignores_noise():
    import pytest
    from molorient.helpers.diagonalization import rot_mat_from_axes
    from molorient.classes.vector import Vector
    from molorient.orientation.standardize_axes.axis_standardization import pick_canonical_frame

    def unit(i, sign=1):
        v = Vector(3)
        v.assign(i, sign)
        return v

    atoms = [Atom("H", 0, 0, Decimal('0.5'), 1), Atom("H", 0, 0, Decimal('-0.5'), 1),
             Atom("O", 0, 0, Decimal('1e-30'), 8)]
    proper = rot_mat_from_axes([unit(0), unit(1), unit(2)])
    flipped = rot_mat_from_axes([unit(0), unit(1, -1), unit(2, -1)])
    tol = Decimal('1e-6')
    #both senses of the axis give the same atoms within tol; the first is returned
    assert pick_canonical_frame([proper, flipped], atoms, tol) is proper
    assert pick_canonical_frame([flipped, proper], atoms, tol) is flipped
    mirror = rot_mat_from_axes([unit(0), unit(1), unit(2, -1)])
    with pytest.raises(ValueError):
        pick_canonical_frame([mirror], atoms, tol)


def test_cn_axes_finder_raises_without_axes():
    import pytest
    atoms = [Atom("C", 0, 0, 0, 6), Atom("H", 1, 0, 0, 1), Atom("H", 0, 2, 0, 1), Atom("H", 0, 0, 3, 1)]
    with pytest.raises(ValueError):
        cn_axes_finder(atoms)
