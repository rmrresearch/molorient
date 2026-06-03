import numpy as np
from molorient.classes.atom import Atom
from molorient.utils.axis_standardization import inertia_tensor, standardize_axes

def test_axis_standardization():
    # #Test case for inertia tensor.
    # atoms = [
    #     Atom("H", 0.0, 0.0, 1.0, 1.008, 1.0),
    #     Atom("H", 0.0, 1.0, 0.0, 1.008, 1.0),
    #     Atom("H", 1.0, 0.0, 1.0, 1.008, 1.0),
    # ]

    # moment_a, moment_b, moment_c, norm_eigvecs = inertia_tensor(atoms)

    # assert moment_a < moment_b < moment_c
    # assert np.isclose(moment_a, (5-np.sqrt(5))/2)
    # assert np.isclose(moment_b, 3.0)
    # assert np.isclose(moment_c, (5+np.sqrt(5))/2)

    # #Test for a single atom
    # atoms_single = [
    #     Atom("H", 0.0, 0.0, 0.0, 1.008, 1.0),
    # ]

    # moment_a, moment_b, moment_c, norm_eigvecs = inertia_tensor(atoms_single)
    # standardized_atoms = standardize_axes(atoms_single, moment_a, moment_b, moment_c, norm_eigvecs)
    # for atom in standardized_atoms:
    #     assert np.isclose(float(atom.x), 0.0) 
    #     assert np.isclose(float(atom.y), 0.0)
    #     assert np.isclose(float(atom.z), 0.0)

    # #Test case for an asymmetric top molecule
    # atoms_asymm = [
    #     Atom("H", 0.0, 0.0, 1.0, 1.008, 1.0),
    #     Atom("H", 0.0, 1.0, 0.0, 1.008, 1.0),
    #     Atom("H", 1.0, 0.0, 1.0, 1.008, 1.0),
    # ]

    # moment_a, moment_b, moment_c, norm_eigvecs = inertia_tensor(atoms_asymm)
    # standardized_atoms = standardize_axes(atoms_asymm, moment_a, moment_b, moment_c, norm_eigvecs)
    # expected_asymm_coords = [
    #     (0.851, 0.000, -0.526),
    #     (0.000, 1.000, 0.000),
    #     (1.377, 0.000, 0.325),
    # ]

    # for atom, (x_exp, y_exp, z_exp) in zip(standardized_atoms, expected_asymm_coords):
    #     assert np.isclose(float(atom.x), x_exp, atol=1e-3)
    #     assert np.isclose(float(atom.y), y_exp, atol=1e-3)
    #     assert np.isclose(float(atom.z), z_exp, atol=1e-3)

    # #Test case for a linear molecule
    # atoms_linear = [
    #     Atom("H", 0.0, 0.0, 0.0, 1.008, 1.0),
    #     Atom("H", 0.0, 1.0, 0.0, 1.008, 1.0),
    # ]

    # moment_a, moment_b, moment_c, norm_eigvecs = inertia_tensor(atoms_linear)
    # standardized_atoms = standardize_axes(atoms_linear, moment_a, moment_b, moment_c, norm_eigvecs)
    # expected_linear_coords = [
    #     (0.0, 0.0, 0.0),
    #     (0.0, 0.0, -1.0),
    # ]

    # for atom, (x_exp, y_exp, z_exp) in zip(standardized_atoms, expected_linear_coords):
    #     assert np.isclose(float(atom.x), x_exp, atol=1e-3)
    #     assert np.isclose(float(atom.y), y_exp, atol=1e-3)
    #     assert np.isclose(float(atom.z), z_exp, atol=1e-3)

    #Test case for a symmetric top molecule
    atoms_symm = [
        Atom("He", 0.0, 0.0, 0.0, 4.0026, 2.0),
        Atom("H", 0.0, 0.0, 1.0, 1.008, 1.0),
        Atom("H", 0.0, -0.8660254, -0.5, 1.008, 1.0),
        Atom("H", 0.0, 0.8660254, -0.5, 1.008, 1.0),
    ]

    moment_a, moment_b, moment_c, norm_eigvecs = inertia_tensor(atoms_symm)
    standardized_atoms = standardize_axes(atoms_symm, moment_a, moment_b, moment_c, norm_eigvecs)
    expected_symm_coords = [
     (0.0, 0.0, 0.0),
     (0.0, 1.0, 0.0),
     (-0.8660254, -0.5, 0.0),
     (0.8660254, -0.5, 0.0),       
    ]

    for atom, (x_exp, y_exp, z_exp) in zip(standardized_atoms, expected_symm_coords):
        assert np.isclose(float(atom.x), x_exp, atol=1e-3)
        assert np.isclose(float(atom.y), y_exp, atol=1e-3)
        assert np.isclose(float(atom.z), z_exp, atol=1e-3)