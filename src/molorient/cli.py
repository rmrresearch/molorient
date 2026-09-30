from molorient.classes.atom import Atom
from molorient.orient_system import orient_system
from molorient.orientation.standardize_axes.nwchem_orientation import orient_mol
import periodictable as pt
from decimal import getcontext
import argparse
import os


_MODES = {
    "molorient": (orient_system, "Standardized geometry by molorient"),
    "nwchem": (orient_mol, "Standardized geometry by molorient (NWChem orientation)"),
}


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


def write_xyz(filepath, atoms, comment):
    """Writes atoms to filepath in .xyz format, with comment as the second line."""
    with open(filepath, 'w') as f:
        f.write(f"{len(atoms)}\n")
        f.write(f"{comment}\n")
        for atom in atoms:
            f.write(f"{atom.element} {atom.x} {atom.y} {atom.z}\n")


def main():
    """
    Runs parse_xyz(), set_precision(), and orient_system() (or orient_mol(),
    for NWChem orientation) to standardize geometry.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("file", help="Path to file")
    parser.add_argument(
        "mode", nargs="?", default="molorient", type=str.lower,
        choices=list(_MODES),
        help="Orientation to use: 'molorient' (default) or 'nwchem'"
    )
    args = parser.parse_args()
    atoms, folder, base, ext = parse_xyz(args.file)
    set_precision()

    orient, comment = _MODES[args.mode]
    std_atoms = orient(atoms)

    std_xyz_filepath = os.path.join(folder, f"{base}_standardized{ext}")
    write_xyz(std_xyz_filepath, std_atoms, comment)
