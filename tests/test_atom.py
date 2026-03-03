from atom import Atom
from decimal import Decimal

def test_atom_initialization():
    atom = Atom(
        "C",
        Decimal("1.0"),
        Decimal("2.0"),
        Decimal("3.0"),
        Decimal("12.011"),
        Decimal("6")
    )
    
    assert atom.element == "C"
    assert atom.x == Decimal("1.0")
    assert atom.y == Decimal("2.0")
    assert atom.z == Decimal("3.0")
    assert atom.mass == Decimal("12.011")
    assert atom.charge == Decimal("6")