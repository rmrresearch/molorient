from molorient.utils.orient_system import orient_system
from molorient.utils.cli import parse_xyz, set_precision
import argparse
import os


def main():
    """
    Runs parse_xyz(), set_precision(), and orient_system() to standardize geometry.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("file", help="Path to file")
    args = parser.parse_args()
    atoms, folder, base, ext = parse_xyz(args.file)
    set_precision()
    std_atoms = orient_system(atoms)
    std_xyz_filepath = os.path.join(folder, f"{base}_standardized{ext}")

    with open(std_xyz_filepath, 'w') as f:
        f.write(f"{len(std_atoms)}\n")
        f.write("Standardized geometry by molorient\n")
        for atom in std_atoms:
            f.write(f"{atom.element} {atom.x} {atom.y} {atom.z}\n")


if __name__ == '__main__':
    main()