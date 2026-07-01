from molorient.utils.randomize_utils.random_permutations import random_permutation
from molorient.utils.cli import parse_xyz, set_precision
from molorient.utils.orient_system import orient_system
from unittest.mock import patch


def test_permute_asymmetric(tmp_path):
    test_file = tmp_path / "test_asymm.xyz"
    test_file.write_text(
    "11\n"
    "Comment line\n"
    "H -2.16595 0.000000 0.355307\n"
    "H -1.30033 -0.881697 -0.905080\n"
    "H -1.30033 0.881697 -0.905080\n"
    "H 0.000000 -0.875477 1.24649\n"
    "H 0.000000 0.875477 1.24649\n"
    "H 1.30033 -0.881697 -0.905080\n"
    "H 1.30033 0.881697 -0.905080\n"
    "H 2.16595 0.000000 0.355307\n"
    "C -1.26529 0.000000 -0.260762\n"
    "C 0.000000 0.000000 0.590977\n"
    "C 1.26529 0.000000 -0.260762\n"
    )
    
    atoms, folder, base, ext = parse_xyz(test_file)
    with patch('builtins.input', return_value = '6'):
        result = set_precision()
        assert result == 6
    shuffled_atoms = random_permutation(atoms)
    std_atoms = orient_system(shuffled_atoms)

    for a in sorted(atoms, key=lambda a: (a.element, a.x, a.y, a.z)):
        print("REF", a.element, a.x, a.y, a.z)
    for a in sorted(std_atoms, key=lambda a: (a.element, a.x, a.y, a.z)):
        print("STD", a.element, a.x, a.y, a.z)

    for a in std_atoms:
        print(a.element, a.x, a.y, a.z)

    for atom, ref in zip(
        sorted(std_atoms, key=lambda a: (a.element, a.x, a.y, a.z)),
        sorted(atoms, key=lambda a: (a.element, a.x, a.y, a.z)),
    ):
        assert atom.element == ref.element
        assert atom.x == ref.x
        assert atom.y == ref.y
        assert atom.z == ref.z