from decimal import Decimal, getcontext
from molorient.utils.axis_standardization import fix_molecule_sign


def sort_atoms(atoms, eigvals):
    """
    Sorts standardized atoms. Lighter elements go first, followed by lowest x coordinate, then lowest y coordinate, 
    then lowest z coordinate.
    """

    # Round atom coordinates
    # (1) If number of decimal places > number of sig figs, round to 10**(-precision)
    # (2) If number of sig figs > number of decimal places, round to nearest precision
    for atom in atoms:
        coords = [atom.x, atom.y, atom.z]

        for i, coord in enumerate(coords):
            t = coord.as_tuple()
            sig_figs = len(t.digits)
            dec_places = max(0, -t.exponent)

            if sig_figs < dec_places:
                coords[i] = coord.quantize(Decimal(10) ** -(getcontext().prec))
            else:
                rounded = round(coord, getcontext().prec - coord.adjusted() - 1)
                coords[i] = Decimal(str(rounded))

        atom.x = coords[0]
        atom.y = coords[1]
        atom.z = coords[2]

    atoms.sort(key = lambda atom: (atom.charge, atom.x, atom.y, atom.z))

    #Deals with asymmetric top case:
    if eigvals[0] != eigvals[1] != eigvals[2]:
        atoms = fix_molecule_sign(atoms)
        atoms.sort(key = lambda atom: (atom.charge, atom.x, atom.y, atom.z))

    return atoms