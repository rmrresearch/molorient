from molorient.utils.orient_system import orient_system
from molorient.classes.atom import Atom
from decimal import Decimal, getcontext

def test_orient_system():
    getcontext().prec = 28
    atoms = [
        Atom("H", 0, 0, 0, 1),
        Atom("H", 1, 1, 1, 1),
        Atom("H", 1, -1, 1, 1)
    ]

    std_atoms = orient_system(atoms)
    assert std_atoms[0].x == -1 
    assert std_atoms[0].y == Decimal('-0.4714045207910316829338962414')
    assert std_atoms[0].z == 0

    assert std_atoms[1].x == 0
    assert std_atoms[1].y == Decimal('0.9428090415820633658677924828')
    assert std_atoms[1].z == 0

    assert std_atoms[2].x == 1
    assert std_atoms[2].y == Decimal('-0.4714045207910316829338962414')
    assert std_atoms[2].z == 0