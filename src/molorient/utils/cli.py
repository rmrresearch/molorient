from molorient.classes.atom import Atom
import periodictable as pt
from decimal import getcontext
import os


def parse_xyz(filepath):
    """
    Parses the .xyz file for element and coordinates and assigns charge.
    """
    atoms = []

    with open(filepath, 'r') as f:
        lines = f.readlines()

    for l in lines[2:]:
        l = l.strip()
        if not l:
            continue
        parts = l.split()
        element = parts[0]
        x, y, z = parts[1], parts[2], parts[3]
        el = pt.elements.symbol(element)
        #Assigns charge
        charge = el.number
        atoms.append(Atom(element, x, y, z, charge))

    folder = os.path.dirname(filepath)
    filename = os.path.basename(filepath)
    base, ext = os.path.splitext(filename)

    return atoms, folder, base, ext


def set_precision():
    """
    Allows the user to specify number of sig figs for calculations.
    """
    user_prec = int(input("Enter desired number of significant figures: "))
    getcontext().prec = user_prec

    return user_prec