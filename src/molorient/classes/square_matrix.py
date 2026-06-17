from decimal import Decimal
from molorient.classes.vector import Vector
from  decimal import Decimal, getcontext


class SquareMatrix:
    def __init__(self, n):
        self.elements = [[Decimal('0')] * n for _ in range(n)]
    
    def assign(self, i, j, value):
        self.elements[i][j] = Decimal(value)
    
    def add(self, other):
        n = len(self.elements)
        result = SquareMatrix(n)
        for i in range(n):
            for j in range(n):
                result.elements[i][j] = self.elements[i][j] + other.elements[i][j]
        return result

    def multiply(self, other):
        n = len(self.elements)
        if isinstance(other, SquareMatrix):
            result = SquareMatrix(n)
            for i in range(n):
                for j in range(n):
                    for k in range(n):
                        result.elements[i][j] += self.elements[i][k] * other.elements[k][j]
            return result
    
        elif isinstance(other, Vector):
            result = Vector(n)
            for i in range(n):
                    for j in range(n):
                        result.elements[i] += self.elements[i][j] * other.elements[j]
            return result
        
    def transpose(self):
        n = len(self.elements)
        result = SquareMatrix(n)
        for i in range(n):
            for j in range(n):
                result.elements[i][j] = self.elements[j][i]
        return result
    
    def scale(self, scalar):
        n = len(self.elements)
        result = SquareMatrix(n)
        for i in range(n):
            for j in range(n):
                result.elements[i][j] = Decimal(scalar) * self.elements[i][j]
        return result