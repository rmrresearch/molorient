from decimal import Decimal


class Atom:
    def __init__(self, element, x, y, z, charge):
        self.element = element
        self.x = Decimal(x)
        self.y = Decimal(y)
        self.z = Decimal(z)
        self.charge = Decimal(charge)