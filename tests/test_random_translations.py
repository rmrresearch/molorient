from molorient.utils.cli import parse_xyz, set_precision
from molorient.utils.orient_system import orient_system
from unittest.mock import patch
from decimal import Decimal, getcontext
import random


def test_translate_asymmetric():
    """
    Randomly translates urea_ref.xyz and checks that randomly translating atoms gives same output geometry.
    """

    atoms, folder, base, ext = parse_xyz("tests/randomized_tests/translations/urea_ref.xyz")
    oriented_ref_atoms, folder, base, ext = parse_xyz("tests/randomized_tests/translations/urea_oriented_6.xyz")
    with patch('builtins.input', return_value = '6'):
        result = set_precision()
        assert result == 6
    getcontext().prec += 5

    trans_vec = [Decimal(random.randint(-1000000,1000000))/Decimal(100000) for _ in range(3)]    
    translated_atoms = atoms

    for atom in translated_atoms:
        atom.x += trans_vec[0]
        atom.y += trans_vec[1]
        atom.z += trans_vec[2]

    getcontext().prec -= 5
    std_atoms = orient_system(translated_atoms)

    for atom, ref in zip(std_atoms, oriented_ref_atoms):
        assert atom.element == ref.element
        assert atom.x == ref.x
        assert atom.y == ref.y
        assert atom.z == ref.z


def test_translate_symmetric():
    """
    Randomly translates benzene_ref.xyz and checks that randomly translating atoms gives same output geometry.
    """

    atoms, folder, base, ext = parse_xyz("tests/randomized_tests/translations/benzene_ref.xyz")
    oriented_ref_atoms, folder, base, ext = parse_xyz("tests/randomized_tests/translations/benzene_oriented_6.xyz")

    with patch('builtins.input', return_value = '6'):
        result = set_precision()
        assert result == 6

    for atom in atoms:
        atom.x = +atom.x
        atom.y = +atom.y
        atom.z = +atom.z

    getcontext().prec += 5

    trans_vec = [Decimal(random.randint(-1000000,1000000))/Decimal(100000) for _ in range(3)]
    translated_atoms = atoms

    for atom in translated_atoms:
        atom.x += trans_vec[0]
        atom.y += trans_vec[1]
        atom.z += trans_vec[2]

    getcontext().prec -= 5
    std_atoms = orient_system(translated_atoms)

    for atom, ref in zip(std_atoms, oriented_ref_atoms):
        assert atom.element == ref.element
        assert atom.x == ref.x
        assert atom.y == ref.y
        assert atom.z == ref.z


def test_translate_spherical():
    """
    Randomly translates sf6_ref.xyz and checks that randomly translating atoms gives same output geometry.
    """

    atoms, folder, base, ext = parse_xyz("tests/randomized_tests/translations/sf6_ref.xyz")
    oriented_ref_atoms, folder, base, ext = parse_xyz("tests/randomized_tests/translations/sf6_oriented_6.xyz")

    with patch('builtins.input', return_value = '6'):
        result = set_precision()
        assert result == 6

    for atom in atoms :
        atom.x = +atom.x
        atom.y = +atom.y
        atom.z = +atom.z   
    
    getcontext().prec += 5
    trans_vec = [Decimal(random.randint(-1000000,1000000))/Decimal(100000) for _ in range(3)]   
    translated_atoms = atoms

    for atom in translated_atoms:
        atom.x += trans_vec[0]
        atom.y += trans_vec[1]
        atom.z += trans_vec[2]

    getcontext().prec -= 5
    std_atoms = orient_system(translated_atoms)

    for atom, ref in zip(std_atoms, oriented_ref_atoms):
        assert atom.element == ref.element
        assert atom.x == ref.x
        assert atom.y == ref.y
        assert atom.z == ref.z