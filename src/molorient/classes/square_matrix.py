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
    
    def row_reduce(self):
        tol = 10 ** -(Decimal(getcontext().prec - 2))

        n = len(self.elements)
        result = SquareMatrix(n)
        for i in range(n):
            for j in range(n):
                result.elements[i][j] = self.elements[i][j]

        pivot_row = 0
        for col in range(n):
            max_val = abs(result.elements[pivot_row][col])
            max_row = pivot_row
            for r in range(pivot_row + 1, n):
                if abs(result.elements[r][col]) > max_val:
                    max_val = abs(result.elements[r][col])
                    max_row = r

            if max_val < tol:
                continue

            if max_row != pivot_row:
                result.elements[pivot_row], result.elements[max_row] = (
                    result.elements[max_row], result.elements[pivot_row]
                )
            
            pivot_val = result.elements[pivot_row][col]
            for c in range(n):
                result.elements[pivot_row][c] = result.elements[pivot_row][c] / pivot_val
            
            for r in range(n):
                if r != pivot_row:
                    factor = result.elements[r][col]
                    if abs(factor) > tol:
                        for c in range(n):
                            if c == col:
                                result.elements[r][c] = Decimal('0')
                            else:
                                result.elements[r][c] -= factor * result.elements[pivot_row][c]
            
            pivot_row += 1
            if pivot_row == n:
                break

        return result