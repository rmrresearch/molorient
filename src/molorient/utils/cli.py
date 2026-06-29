from molorient.classes.atom import Atom
from molorient.utils.orient_system import orient_system
import argparse
import periodictable as pt
from decimal import getcontext
import os


def read_xyz(filepath):
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


def get_precision():
    """
    Allows the user to specify number of sig figs.
    """
    user_prec = int(input("Enter desired number of significant figures: "))
    getcontext().prec = user_prec

    return user_prec


def main():
    """
    Runs read_xyz(), get_precision(), and orient_system() to standardize geometry.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("file", help="Path to file")
    args = parser.parse_args()
    atoms, folder, base, ext = read_xyz(args.file)
    get_precision()
    std_atoms = orient_system(atoms)
    std_xyz_filepath = os.path.join(folder, f"{base}_standardized{ext}")

    with open(std_xyz_filepath, 'w') as f:
        f.write(f"{len(std_atoms)}\n")
        f.write("Standardized geometry by molorient\n")
        for atom in std_atoms:
            f.write(f"{atom.element} {atom.x} {atom.y} {atom.z}\n")


if __name__ == '__main__':
    main()