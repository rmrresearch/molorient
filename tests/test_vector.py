from molorient.classes.vector import Vector
from decimal import Decimal
import numpy as np
n = 3


def test_vector_initialization():
    zero_vec = Vector(n)
    np_vec = np.zeros(n)

    for i in range(n):
        assert zero_vec.elements[i] == Decimal(np_vec[i])


def test_assign():
    vec = Vector(n)

    for i in range(n):
        vec.assign(i, Decimal('1'))

    arr = np.array([1.0, 1.0, 1.0])

    for i in range(n):
        assert vec.elements[i] == Decimal(arr[i])


def test_add():
    u = Vector(n)
    v = Vector(n)

    for i in range(n):
        u.elements[i] = Decimal('4.0')
        v.elements[i] = Decimal('3.0')
    
    w = u.add(v)

    arr = np.array([7.0, 7.0, 7.0])     

    for i in range(n):
        assert  w.elements[i] == Decimal(arr[i])


def test_dot():
    u = Vector(n)
    v = Vector(n)

    for i in range(n):
        u.elements[i] = Decimal('2.0')
        v.elements[i] = Decimal('3.0')
    
    x = u.dot(v)

    assert x == Decimal('18.0')


def test_scale():
    v = Vector(n)

    for i in range(n):
        v.elements[i] = Decimal('1.0')
    
    w = v.scale(Decimal('3.0'))

    arr = np.array([3.0, 3.0, 3.0])

    for i in range(n):
        assert w.elements[i] == Decimal(arr[i])


def test_cross():
    v = Vector(n)
    u = Vector(n)

    for i in range(n):
        v.elements[i] = Decimal('2.0')
        u.elements[i] = Decimal('7.0')
    
    u.assign(1, Decimal('0.0'))
    
    w = v.cross(u)

    arr = np.array([14.0, 0.0, -14.0])
    
    for i in range(n):
        assert w.elements[i] == Decimal(arr[i])

