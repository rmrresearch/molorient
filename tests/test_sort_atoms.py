from molorient.classes.atom import Atom
from molorient.utils.sort_atoms import sort_atoms


def test_sort_atoms():
    atoms = [
        Atom("C", 0, 2, 1, 6),
        Atom("H", 0, 2, 0, 1),
        Atom("H", 1, 2, 3, 1),
        Atom("Sb", 0, 2, 3, 51)
    ]

    sorted_atoms = sort_atoms(atoms)

    assert sorted_atoms[0].element == 'H'
    assert sorted_atoms[0].x == 0
    assert sorted_atoms[1].element == 'H'
    assert sorted_atoms[2].element == 'C'
    assert sorted_atoms[3].element == 'Sb'