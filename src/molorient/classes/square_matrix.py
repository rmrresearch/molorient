from decimal import Decimal
from molorient.classes.vector import Vector


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

        raise TypeError(f"Cannot multiply SquareMatrix by {type(other).__name__}.")

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

    def negate(self):
        n = len(self.elements)
        result = SquareMatrix(n)
        for i in range(n):
            for j in range(n):
                result.elements[i][j] = -self.elements[i][j]
        return result

    @classmethod
    def identity(cls, n):
        result = cls(n)
        for i in range(n):
            result.elements[i][i] = Decimal('1')
        return result

    @classmethod
    def from_columns(cls, columns):
        """columns[0], columns[1], columns[2] become this matrix's X, Y, Z columns."""
        n = len(columns)
        result = cls(n)
        result.elements = [[columns[col].elements[row] for col in range(n)] for row in range(n)]
        return result

    def inverse(self):
        """
        Only correct when self is orthogonal (columns orthonormal), where
        A^-1 == A^T -- e.g. an eigenvector matrix. Not a general 3x3 inverse.
        """
        n = len(self.elements)
        result = SquareMatrix(n)

        a = Vector(n)
        b = Vector(n)
        c = Vector(n)

        for i in range(n):
            a.elements[i] = self.elements[i][0]
            b.elements[i] = self.elements[i][1]
            c.elements[i] = self.elements[i][2]

        det = a.dot(b.cross(c))
        transpose = self.transpose()

        result = transpose.scale(Decimal('1') / det)

        return result