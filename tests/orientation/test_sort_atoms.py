from molorient.classes.atom import Atom
from molorient.utils.sort_atoms import sort_atoms
from molorient.utils.axis_standardization import inertia_tensor

def test_sort_atoms():
    atoms = [
        Atom("C", 0, 2, 1, 6),
        Atom("H", 0, 2, 0, 1),
        Atom("H", 1, 2, 3, 1),
        Atom("Sb", 0, 2, 3, 51)
    ]
    moments, _ = inertia_tensor(atoms)

    sorted_atoms = sort_atoms(atoms, moments)

    assert sorted_atoms[0].element == 'H'
    assert sorted_atoms[0].x == -1
    assert sorted_atoms[1].element == 'H'
    assert sorted_atoms[2].element == 'C'
    assert sorted_atoms[3].element == 'Sb'