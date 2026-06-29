from molorient.classes.atom import Atom
import argparse
import periodictable as pt



def read_xyz(filepath):

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
        charge = el.number
        atoms.append(Atom(element, x, y, z, charge))

    return atoms


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("file", help="Path to file")
    args = parser.parse_args()
    return read_xyz(args.file)