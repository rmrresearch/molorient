from decimal import Decimal
from molorient.classes.vector import Vector


class Atom:
    def __init__(self, element, x, y, z, charge):
        self.element = element
        self.x = Decimal(x)
        self.y = Decimal(y)
        self.z = Decimal(z)
        self.charge = Decimal(charge)

    def __repr__(self):
        return f"Atom({self.element!r}, {self.x}, {self.y}, {self.z}, {self.charge})"

    def __eq__(self, other):
        if not isinstance(other, Atom):
            return NotImplemented
        return (self.element, self.x, self.y, self.z, self.charge) == \
               (other.element, other.x, other.y, other.z, other.charge)

    def __hash__(self):
        return hash((self.element, self.x, self.y, self.z, self.charge))

    def position(self):
        """This atom's (x, y, z) as a Vector(3)."""
        return Vector.from_atom(self)