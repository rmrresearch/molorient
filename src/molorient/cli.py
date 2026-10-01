from molorient.classes.atom import Atom
from molorient.orient_system import orient_system
from molorient.orientation.standardize_axes.nwchem_orientation import orient_mol
import periodictable as pt
from decimal import getcontext
import argparse
import os


_MODES = {
    "molorient": (orient_system, "Standardized geometry by MolOrient"),
    "nwchem": (orient_mol, "Standardized geometry by MolOrient (NWChem orientation)"),
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


def set_precision(prec):
    """
    Sets the number of sig figs used for calculations.
    """
    getcontext().prec = prec

    return prec


def positive_int(value):
    """
    Argparse type for --precision: an integer greater than zero.
    """
    try:
        number = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"invalid integer: {value!r}")
    if number <= 0:
        raise argparse.ArgumentTypeError("precision must be greater than 0")

    return number


def build_parser():
    """
    Builds the command-line parser.
    """
    parser = argparse.ArgumentParser(prog="molorient")
    parser.add_argument("file", nargs="?", default=None,
                        help="Path to .xyz file (same as --input)")
    parser.add_argument("--input", dest="input_file", default=None,
                        help="Path to .xyz file")
    orientation = parser.add_mutually_exclusive_group()
    orientation.add_argument("--sno", action="store_true",
                             help="Standard Nuclear Orientation (default)")
    orientation.add_argument("--nwchem", action="store_true",
                             help="NWChem orientation")
    parser.add_argument("--precision", type=positive_int, default=6,
                        help="Number of significant figures (default: 6)")

    return parser


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
    for --nwchem) to standardize geometry.
    """
    parser = build_parser()
    args = parser.parse_args()

    if args.file is None and args.input_file is None:
        parser.error("an input .xyz file is required (positional or --input)")
    if (args.file is not None and args.input_file is not None
            and args.file != args.input_file):
        parser.error("positional file and --input disagree; give only one")
    filepath = args.input_file if args.input_file is not None else args.file

    atoms, folder, base, ext = parse_xyz(filepath)
    set_precision(args.precision)

    orient, comment = _MODES["nwchem" if args.nwchem else "molorient"]
    std_atoms = orient(atoms)

    std_xyz_filepath = os.path.join(folder, f"{base}_standardized{ext}")
    write_xyz(std_xyz_filepath, std_atoms, comment)
