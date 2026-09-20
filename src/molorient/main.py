from molorient.utils.orient_system import orient_system
from molorient.utils.nwchem_orientation import orient_mol
from molorient.utils.cli import parse_xyz, set_precision, write_xyz
import argparse
import os


_MODES = {
    "molorient": (orient_system, "Standardized geometry by molorient"),
    "nwchem": (orient_mol, "Standardized geometry by molorient (NWChem orientation)"),
}


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


if __name__ == '__main__':
    main()
