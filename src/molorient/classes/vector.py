from decimal import Decimal


class Vector:
    def __init__(self, n):
        self.elements = [Decimal('0')] * n
    
    def assign(self, i, value):
        self.elements[i] = Decimal(value)
    
    def add(self, other):
        n = len(self.elements)
        result = Vector(n)
        for i in range(n):
            result.elements[i] = self.elements[i] + other.elements[i]
        return result

    def subtract(self, other):
        n = len(self.elements)
        result = Vector(n)
        for i in range(n):
            result.elements[i] = self.elements[i] - other.elements[i]
        return result

    def dot(self, other):
        if len(self.elements) != len(other.elements):
            raise ValueError("Vectors must be the same dimension.")
        n = len(self.elements)
        result = Decimal('0')
        for i in range(n):
            result += self.elements[i] * other.elements[i]
        return result

    def scale(self, scalar):
        n = len(self.elements)
        scalar = Decimal(scalar)
        result = Vector(n)
        for i in range(n):
            result.elements[i] = scalar * self.elements[i]
        return result

    def negate(self):
        n = len(self.elements)
        result = Vector(n)
        for i in range(n):
            result.elements[i] = -self.elements[i]
        return result

    def norm(self):
        return self.dot(self).sqrt()

    def normalized(self):
        return self.scale(1 / self.norm())

    @classmethod
    def from_atom(cls, atom):
        """A Vector(3) holding an atom's (x, y, z)."""
        v = cls(3)
        v.assign(0, atom.x)
        v.assign(1, atom.y)
        v.assign(2, atom.z)
        return v

    def cross(self, other):
        if len(self.elements) != 3 or len(other.elements) != 3:
            raise ValueError("cross() is only defined for 3-dimensional vectors.")
        result = Vector(3)
        result.elements[0] = (self.elements[1] * other.elements[2]) - (self.elements[2] * other.elements[1])
        result.elements[1] = (self.elements[2] * other.elements[0]) - (self.elements[0] * other.elements[2])
        result.elements[2] = (self.elements[0] * other.elements[1]) - (self.elements[1] * other.elements[0])
        return result