from molorient.classes.square_matrix import SquareMatrix
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


def test_multiply():
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
