from molorient.utils.cli import parse_xyz, set_precision
from molorient.utils.orient_system import orient_system
from molorient.utils.axis_standardization import inertia_tensor
from unittest.mock import patch
import random



def test_permute_asymmetric():
    """
    Randomly permutes butane_ref.xyz and checks that randomly permutating atoms
    gives same output geometry.
    """
    atoms, folder, base, ext = parse_xyz("tests/randomized_tests/permutations/butane_ref.xyz")
    oriented_ref_atoms, folder, base, ext = parse_xyz("tests/randomized_tests/permutations/butane_oriented_6.xyz")
    with patch('builtins.input', return_value = '6'):
        result = set_precision()
        assert result == 6
    shuffled_atoms = random.sample(atoms, len(atoms))
    std_atoms = orient_system(shuffled_atoms)

    for atom, ref in zip(std_atoms, oriented_ref_atoms):
        assert atom.element == ref.element
        assert atom.x == ref.x
        assert atom.y == ref.y
        assert atom.z == ref.z