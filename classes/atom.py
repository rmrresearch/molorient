from decimal import Decimal

class Atom:
    def __init__(self, element, x, y, z, mass, charge):
        self.element = element
        self.x = Decimal(x)
        self.y = Decimal(y)
        self.z = Decimal(z)
        self.mass = Decimal(mass)
        self.charge = Decimal(charge)