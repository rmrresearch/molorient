def sort_atoms(atoms):
    """
    Sorts standardized atoms. Lighter elements go first, followed by lowest x coordinate, then lowest y coordinate, 
    then lowest z coordinate.
    """

    atoms.sort(key = lambda atom: (atom.charge, atom.x, atom.y, atom.z))

    for atom in atoms:
        atom.x = +atom.x
        atom.y = +atom.y
        atom.z = +atom.z

    return atoms