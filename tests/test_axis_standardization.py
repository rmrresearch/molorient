from molorient.utils.axis_standardization import inertia_tensor
from molorient.utils.diagonalization import eigval_solver
from molorient.classes.atom import Atom
import numpy as np
from decimal import Decimal


def test_inertia_tensor():
    atoms =[
        Atom("H", 0.0, 0.0, 0.0, 1.008, 1.0),
        Atom("H", 1.0, 0.0, 0.0, 1.008, 1.0)
    ]

    moments, eigvecs = inertia_tensor(atoms) 

    vecs = np.array([
        [1.0, 0.0, 0.0],
        [0.0, 0.0, -1.0],
        [0.0, 1.0 ,0.0]
    ])

    tol = Decimal('10')**(Decimal('-6'))
    assert (moments[0] - 0) < tol
    assert (moments[1] - 1) < tol
    assert (moments[2] - 1) < tol
    for i in range(3):
        assert eigvecs[0].elements[i] - Decimal(vecs[0][i]) < tol
        assert eigvecs[1].elements[i] - Decimal(vecs[1][i]) < tol
        assert eigvecs[2].elements[i] - Decimal(vecs[2][i]) < tol