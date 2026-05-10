import numpy as np
from molorient.classes.atom import Atom
from molorient.utils.axis_standardization import inertia_tensor, standardize_axes

def test_axis_standardization():
    #Test case for inertia tensor.
    atoms = [
        Atom("H", 0.0, 0.0, 1.0, 1.008, 1.0),
        Atom("H", 0.0, 1.0, 0.0, 1.008, 1.0),
        Atom("H", 1.0, 0.0, 1.0, 1.008, 1.0),
    ]

    moment_a, moment_b, moment_c, norm_eigvecs = inertia_tensor(atoms)

    assert moment_a < moment_b < moment_c
    assert np.isclose(moment_a, (5-np.sqrt(5))/2)
    assert np.isclose(moment_b, 3.0)
    assert np.isclose(moment_c, (5+np.sqrt(5))/2)


    # #Test for a single atom
    # atoms_single = [
    #     Atom("H", 0.0, 0.0, 0.0, 1.008, 1.0),
    # ]

    # moment_a, moment_b, moment_c, norm_eigvecs = inertia_tensor(atoms_single)
    # standardized_atoms = standardize_axes(atoms_single, moment_a, moment_b, moment_c, norm_eigvecs)
    # for atom in standardized_atoms:
    #     assert atom.x == 0.0
    #     assert atom.y == 0.0
    #     assert atom.z == 0.0 


    # #Test case for an asymmetric top molecule
    # atoms_asymm = [
    #     Atom("H", -0.4, -0.4, -0.5, 1.008, 1.0),
    #     Atom("H", -0.4, -0.4, 0.5, 1.008, 1.0),
    #     Atom("O", 0.1, 0.1, 0.0, 15.999, 1.0),
    # ]

    # expected_asymm_moments = [0.5, 0.8, 1.3]
    


    # moment_a, moment_b, moment_c = inertia_tensor(atoms_asymm)
    # standardized_atoms = standardize_axes(atoms_asymm, moment_a, moment_b, moment_c)
    # for atom, (x_exp, y_exp, z_exp) in zip(standardized_atoms, expected_asymm_coords):
    #     assert atom.x == x_exp
    #     assert atom.y == y_exp
    #     assert atom.z == z_exp

    # #Test case for a single atom
    # atoms_single = [
    #     Atom("H", 0.0, 0.0, 0.0, 1.008, 1.0),
    # ]

    # expected_single_moments = [0.0, 0.0, 0.0]
    # expected_single_coords = [(0.0, 0.0, 0.0)]

    # moment_a, moment_b, moment_c = inertia_tensor(atoms_single)
    # standardized_atoms = standardize_axes(atoms_single, moment_a, moment_b, moment_c)
    # for atom, (x_exp, y_exp, z_exp) in zip(standardized_atoms, expected_single_coords):
    #     assert atom.x == x_exp
    #     assert atom.y == y_exp
    #     assert atom.z == z_exp

    # #Test case for a linear molecule
    # atoms_linear = [
    #     Atom("H", 0.5, 0.0, 0.0, 1.008, 1.0),
    #     Atom("H", -0.5, 0.0, 0.0, 1.008, 1.0),
    # ]

    # expected_linear_moments = [0.0, 0.5, 0.5]
    # expected_linear_coords = [(0.0, 0.0, 0.5), (0.0, 0.0, -0.5)]

    # moment_a, moment_b, moment_c = inertia_tensor(atoms_linear)
    # standardized_atoms = standardize_axes(atoms_linear, moment_a, moment_b, moment_c)
    # for atom, (x_exp, y_exp, z_exp) in zip(standardized_atoms, expected_linear_coords):
    #     assert atom.x == x_exp
    #     assert atom.y == y_exp
    #     assert atom.z == z_exp
