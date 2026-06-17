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
    
    def dot(self, other):
        n = len(self.elements)
        result = Decimal('0')
        for i in range(n):
            result += self.elements[i] * other.elements[i]
        if len(self.elements) != len(other.elements):
            raise ValueError("Vectors must be the same dimension.")
        return result
        
    def scale(self, scalar):
        n = len(self.elements)
        result = Vector(n)
        for i  in range(n):
            result.elements[i] = Decimal(scalar) * self.elements[i]
        return result
    
    def cross(self, other):
        n = len(self.elements)
        result = Vector(n)
        result.elements[0] = (self.elements[1] * other.elements[2]) - (self.elements[2] * other.elements[1])
        result.elements[1] = (self.elements[2] * other.elements[0] - (self.elements[0] * other.elements[2]))
        result.elements[2] = (self.elements[0] * other.elements[1] - (self.elements[1] * other.elements[2]))
        return result