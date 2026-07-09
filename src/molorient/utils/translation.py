from molorient.classes.atom import Atom
from decimal import Decimal, getcontext 
from molorient.classes.vector import Vector


def translation_vector(atoms):
    """
    This function first calculates the center of positive charge (coc) of the system.
    Then, it creates a translation vector that will move the atoms of the system so that
    the coc is at the origin (0, 0, 0).
    """

    getcontext().prec += 5

    #Calculate the center of positive charge (coc)
    total_charge = sum(atom.charge for atom in atoms)
    
    x_center = sum(atom.x * atom.charge for atom in atoms) / total_charge
    y_center = sum(atom.y * atom.charge for atom in atoms) / total_charge
    z_center = sum(atom.z * atom.charge for atom in atoms) / total_charge
    
    coc = Vector(3)
    coc.assign(0, x_center)
    coc.assign(1, y_center)
    coc.assign(2, z_center)

    #Create the translation vector to move atoms so coc is at origin.
    trans_vec = coc.scale(-1)

    getcontext().prec -= 5

    return trans_vec


def translate_to_origin(atoms, trans_vec):
    """
    This function takes a list of atoms and a translation vector, and translates the coordinates of each atom
    by adding the translation vector to the atom's coordinates. It returns a new list of Atom objects with the translated coordinates.
    """
    getcontext().prec += 5

    translated_atoms = []
    for atom in atoms:
        new_x = atom.x + trans_vec.elements[0]
        new_y = atom.y + trans_vec.elements[1]
        new_z = atom.z + trans_vec.elements[2]

        translated_atoms.append(Atom(atom.element, new_x, new_y, new_z, atom.charge))

    getcontext().prec -= 5

    for atom in translated_atoms:
        atom.x = +atom.x
        atom.y = +atom.y
        atom.z = +atom.z

    return translated_atoms