from molorient.utils.diagonalization import eigval_solver, eigvec_solver
from molorient.classes.square_matrix import SquareMatrix
from decimal import Decimal, getcontext, ROUND_HALF_UP
import numpy as np
rng = np.random.default_rng()


def test_diagonalization():
    eigvals = rng.uniform(low = 0, high = 10000000, size = 3)
    e_0, e_1, e_2 = sorted(eigvals)

    d = np.array([[e_0,0, 0 ],
                [0, e_1, 0],
                [0, 0, e_2]])

    a = rng.uniform(low = 0, high = 10000000, size = (3, 3))

    q, r = np.linalg.qr(a)

    b = q @ d @ q.T
    dec_b = SquareMatrix(3)

    for i in range(3):
        for j in range(3):
            dec_b.elements[i][j] = Decimal(str(b[i][j]))

    eigval_0, eigval_1, eigval_2 = sorted(eigval_solver(dec_b))
    tol = Decimal('10')**-Decimal('6')

    assert abs(Decimal(e_0) - eigval_0).quantize(tol, rounding = ROUND_HALF_UP) < tol
    assert abs(Decimal(e_1) - eigval_1).quantize(tol, rounding = ROUND_HALF_UP) < tol
    assert abs(Decimal(e_2) - eigval_2).quantize(tol, rounding = ROUND_HALF_UP) < tol
