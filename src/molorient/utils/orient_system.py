from molorient.utils.translation import translation_vector, translate_to_origin
from molorient.utils.axis_standardization import inertia_tensor, standardize_axes
from molorient.utils.sort_atoms import sort_atoms


def orient_system(atoms):
    """
    Calls translation and axis standardization functions to orient the molecule.
    """
    trans_vec = translation_vector(atoms)
    trans_atoms = translate_to_origin(atoms, trans_vec)
    moments, eigvecs = inertia_tensor(trans_atoms)
    rotated_atoms = standardize_axes(moments, eigvecs, trans_atoms)
    sorted_atoms = sort_atoms(rotated_atoms)

    return sorted_atoms