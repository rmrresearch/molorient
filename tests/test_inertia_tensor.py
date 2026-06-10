from molorient.classes.atom import Atom
from molorient.utils.axis_standardization import inertia_tensor
from decimal import Decimal, getcontext, ROUND_HALF_UP
import numpy as np


def test_inertia_tensor():
    atoms = [
        Atom("H", 0.0, 0.0, 0.0, 1.008, 1.0),
        Atom("H", 1.0, 0.0, 0.0, 1.008, 1.0),
    ]

    moment_a, moment_b, moment_c = inertia_tensor(atoms)

    assert moment_a < moment_b == moment_c
    assert moment_a == 0
    assert moment_b, moment_c == 1