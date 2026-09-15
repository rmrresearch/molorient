from molorient.utils.orient_system import orient_system
from molorient.utils.nwchem_orientation import orient_mol
from molorient.utils.cli import parse_xyz, set_precision
import argparse
import os


def main():
    """
    Runs parse_xyz(), set_precision(), and orient_system() (or orient_mol(),
    for NWChem orientation) to standardize geometry.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("file", help="Path to file")
    parser.add_argument(
        "mode", nargs="?", default="molorient", type=str.lower,
        choices=["molorient", "nwchem"],
        help="Orientation to use: 'molorient' (default) or 'nwchem'"
    )
    args = parser.parse_args()
    atoms, folder, base, ext = parse_xyz(args.file)
    set_precision()

    if args.mode == "nwchem":
        std_atoms = orient_mol(atoms)
        comment = "Standardized geometry by molorient (NWChem orientation)"
    else:
        std_atoms = orient_system(atoms)
        comment = "Standardized geometry by molorient"

    std_xyz_filepath = os.path.join(folder, f"{base}_standardized{ext}")

    with open(std_xyz_filepath, 'w') as f:
        f.write(f"{len(std_atoms)}\n")
        f.write(f"{comment}\n")
        for atom in std_atoms:
            f.write(f"{atom.element} {atom.x} {atom.y} {atom.z}\n")


if __name__ == '__main__':
    main()