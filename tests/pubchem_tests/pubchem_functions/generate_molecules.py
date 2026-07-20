from molorient.utils.trig_helpers import cos_series, sin_series, pi_as_decimal
from molorient.classes.atom import Atom
from molorient.classes.square_matrix import SquareMatrix
from molorient.classes.vector import Vector
from decimal import Decimal, getcontext
import pubchempy as pcp
import random
import os
file_dir = os.path.dirname(__file__)
pubchem_dir = os.path.abspath(os.path.join(file_dir, '..'))


def search_pubchem():
    """
    Randomly selects a compound from the PubChem database using
    PubChemPy. Returns a list of Atom objects and the CID used.
    """
    for _ in range(1000):
        rand_int = random.randint(1, 119000000)
        try:
            c = pcp.Compound.from_cid(rand_int, record_type="3d")
        except (pcp.BadRequestError, pcp.NotFoundError, pcp.PubChemHTTPError):
            continue

        atoms_output = c.to_dict(properties=["atoms"])

        return atoms_output, rand_int
    

def randomize_orientation(atoms):
    """
    Randomly permutes, translates, and rotates a molecule.
    """

    #Randomly permute molecule
    shuffled_atoms = random.sample(atoms, len(atoms))

    #Randomly translate molecule
    translated_atoms = shuffled_atoms
    for atom in translated_atoms:
        atom.x = +atom.x
        atom.y = +atom.y
        atom.z = +atom.z
    getcontext().prec += 5

    trans_vec = [Decimal(random.randint(-100000,100000))/Decimal(10000) for _ in range(3)]    

    for atom in translated_atoms:
        atom.x += trans_vec[0]
        atom.y += trans_vec[1]
        atom.z += trans_vec[2]
    
    getcontext().prec -= 5
    
    #Rotate molecule about randomly generated Tait-Bryan angles
    getcontext().prec += 10
    rot_mat = SquareMatrix(3)
    alpha = Decimal(str(random.random())) * (Decimal(2) * pi_as_decimal())
    beta = Decimal(str(random.random())) * pi_as_decimal()
    gamma = Decimal(str(random.random())) * (Decimal(2) * pi_as_decimal())

    ca = cos_series(alpha)
    sa = sin_series(alpha)

    cb = cos_series(beta)
    sb = sin_series(beta)

    cg = cos_series(gamma)
    sg = sin_series(gamma)

    rot_mat.assign(0, 0, cb * cg)
    rot_mat.assign(0, 1, -cb * sg)
    rot_mat.assign(0, 2, sb)

    rot_mat.assign(1, 0, (ca * sg) + (sa * sb * cg))
    rot_mat.assign(1, 1, (ca * cg) - (sa * sb * sg))
    rot_mat.assign(1, 2, -sa * cb)

    rot_mat.assign(2, 0, (sa * sg) - (ca * sb * cg))
    rot_mat.assign(2, 1, (sa * cg) + (ca * sb * sg))
    rot_mat.assign(2, 2, ca * cb)

    #Rotation
    rotated_atoms = []

    for atom in translated_atoms:
        pos_vec = Vector(3)
        pos_vec.assign(0, atom.x)
        pos_vec.assign(1, atom.y) 
        pos_vec.assign(2, atom.z)

        new_pos = (rot_mat.transpose()).multiply(pos_vec)

        rotated_atoms.append(Atom(atom.element, 
                                        new_pos.elements[0],
                                        new_pos.elements[1],
                                        new_pos.elements[2],
                                        atom.charge))

    getcontext().prec -= 10

    return rotated_atoms