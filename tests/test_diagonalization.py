from molorient.utils.diagonalization import diagonalize_3_by_3
from decimal import Decimal, getcontext, ROUND_HALF_UP
import numpy as np
rng = np.random.default_rng()

def test_eigval_solver():
    eigvals = rng.uniform(low = 0, high = 10000000, size = 3)
    e_0, e_1, e_2 = sorted(eigvals)

    d = np.array([[e_0,0, 0 ],
                [0, e_1, 0],
                [0, 0, e_2]])

    a = rng.uniform(low = 0, high = 10000000, size = (3, 3))

    q, r = np.linalg.qr(a)

    b = q @ d @ q.T
    b_0 = b[0, :].tolist()
    b_1 = b[1, :].tolist()
    b_2 = b[2, :].tolist()

    decimal_b_0 = [Decimal(x) for x in b_0]
    decimal_b_1 = [Decimal(x) for x in b_1]
    decimal_b_2 = [Decimal(x) for x in b_2]

    eigval_0, eigval_1, eigval_2 = sorted(diagonalize_3_by_3(decimal_b_0, decimal_b_1, decimal_b_2))
    tol = Decimal('10')**-Decimal('6')

    assert abs(Decimal(e_0) - eigval_0).quantize(tol, rounding = ROUND_HALF_UP) < tol
    assert abs(Decimal(e_1) - eigval_1).quantize(tol, rounding = ROUND_HALF_UP) < tol
    assert abs(Decimal(e_2) - eigval_2).quantize(tol, rounding = ROUND_HALF_UP) < tol