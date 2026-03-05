from classes.atom import Atom
import numpy as np


def translation_vector(atoms):
    """
    This function first calculates the center of nuclear charge (coc) of the system.
    Then, it creates a translation vector that will move the atoms of the system so that
    the coc is at the origin (0, 0, 0).
    """

    #Calculate the center of nuclear charge (coc)
    total_charge = sum(atom.charge for atom in atoms)
    if total_charge == 0:
        return np.array([0.0, 0.0, 0.0])
    
    x_center = sum(atom.x * atom.charge for atom in atoms) / total_charge
    y_center = sum(atom.y * atom.charge for atom in atoms) / total_charge
    z_center = sum(atom.z * atom.charge for atom in atoms) / total_charge
    
    coc = np.array([x_center, y_center, z_center])

    #Create the translation vector to move atoms so coc is at origin.
    trans_vec = -coc
    return trans_vec


def translate_to_origin(atoms, trans_vec):
    """
    This function takes a list of atoms and a translation vector, and translates the coordinates of each atom
    by adding the translation vector to the atom's coordinates. It returns a new list of Atom objects with the translated coordinates.
    """

    translated_atoms = []
    for atom in atoms:
        new_x = atom.x + trans_vec[0]
        new_y = atom.y + trans_vec[1]
        new_z = atom.z + trans_vec[2]

        translated_atoms.append(Atom(atom.symbol, new_x, new_y, new_z, atom.mass, atom.charge))

    return translated_atoms