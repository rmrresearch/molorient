from molorient.classes.square_matrix import SquareMatrix
from molorient.classes.vector import Vector
import numpy as np
from decimal import Decimal
n = 3


def test_matrix_initialization():
    zero_matrix = SquareMatrix(n)

    numpy_matrix = np.zeros((3, 3))

    for i in range(n):
        for j in range(n):
            assert zero_matrix.elements[i][j] == Decimal(numpy_matrix[i][j])


def test_assign():
    id_matrix = SquareMatrix(n)
    id_matrix.assign(0, 0, Decimal('1'))
    id_matrix.assign(1, 1, Decimal('1'))
    id_matrix.assign(2, 2, Decimal('1'))

    numpy_id = np.eye(3)

    for i in range(n):
        for j in range(n):
            assert id_matrix.elements[i][j] == Decimal(numpy_id[i][j])


def test_add():
    a = SquareMatrix(n)
    b = SquareMatrix(n)

    for i in range(n):
        for j in range(n):
            a.elements[i][j] = Decimal('2.0')
            b.elements[i][j] = Decimal('3.0')

    c = a.add(b)

    arr = np.array([
        [5.0, 5.0, 5.0],
        [5.0, 5.0, 5.0],
        [5.0, 5.0, 5.0],
    ])

    for i in range(n):
        for j in range(n):
            assert c.elements[i][j] == Decimal(arr[i][j])


def test_multiply_matrix():
    a = SquareMatrix(n)
    b = SquareMatrix(n)

    for i in range(n):
        for j in range(n):
            a.elements[i][j] = Decimal('2.0')
            b.elements[i][j] = Decimal('3.0')
    
    c = a.multiply(b)

    arr = np.array([
        [18.0, 18.0, 18.0],
        [18.0, 18.0, 18.0],
        [18.0, 18.0, 18.0]
    ])

    for i in range(n):
        for j in range(n):
            assert c.elements[i][j] == Decimal(arr[i][j])


def test_multiply_vector():
    a = SquareMatrix(n)
    v = Vector(n)

    for i in range(n):
        for j in range(n):
            a.elements[i][j] = Decimal('2.0')
            v.elements[i] = Decimal('4.0')

    b = a.multiply(v)

    arr = np.array([24.0, 24.0, 24.0])
    
    for i in range(n):
        assert b.elements[i] == Decimal(arr[i])


def test_scale():
    a = SquareMatrix(n)

    for i in range(n):
        for j in range(n):
            a.elements[i][j] = Decimal('5.0')
    
    b = a.scale(Decimal('2'))

    arr = np.array([
        [10.0, 10.0, 10.0],
        [10.0, 10.0, 10.0],
        [10.0, 10.0, 10.0],
    ])

    for i in range(n):
        for j in range(n):
            assert b.elements[i][j] == Decimal(arr[i][j])
