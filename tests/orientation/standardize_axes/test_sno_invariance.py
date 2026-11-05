"""
SNO must give the same geometry for any rotation, translation and atom order
of a molecule. One molecule per fixed bug: linear (HCN), symmetric-top Z sign
(NH3), planar symmetric top (BF3) and Td (CH4); inputs exact and rounded to 7
decimals (rounding is what used to break top classification).
"""
import random
from decimal import Decimal as D, getcontext, ROUND_HALF_UP

from molorient.classes.atom import Atom
from molorient.orient_system import orient_system

CHARGE = {"H": 1, "B": 5, "C": 6, "N": 7, "F": 9}
R3 = D(3).sqrt() / 2
A = D("0.626951")
MOLECULES = {
    "HCN": [("H", 0, 0, "-1.064"), ("C", 0, 0, 0), ("N", 0, 0, "1.156")],
    "NH3": [("N", 0, 0, "0.1"), ("H", "0.94", 0, "-0.27"),
            ("H", "-0.47", D("0.94") * R3, "-0.27"), ("H", "-0.47", -D("0.94") * R3, "-0.27")],
    "BF3": [("B", 0, 0, 0), ("F", "1.3", 0, 0),
            ("F", "-0.65", D("1.3") * R3, 0), ("F", "-0.65", -D("1.3") * R3, 0)],
    "CH4": [("C", 0, 0, 0), ("H", A, A, A), ("H", -A, -A, A), ("H", -A, A, -A), ("H", A, -A, -A)],
}


def moved(molecule, rng, decimals):
    """Random rotation (unit quaternion), translation and shuffle, as Atoms."""
    getcontext().prec = 50
    a, b, c, d = (D(repr(rng.gauss(0, 1))) for _ in range(4))
    n = (a*a + b*b + c*c + d*d).sqrt()
    a, b, c, d = a / n, b / n, c / n, d / n
    rot = [[a*a + b*b - c*c - d*d, 2*(b*c - a*d), 2*(b*d + a*c)],
           [2*(b*c + a*d), a*a - b*b + c*c - d*d, 2*(c*d - a*b)],
           [2*(b*d - a*c), 2*(c*d + a*b), a*a - b*b - c*c + d*d]]
    shift = [D(repr(rng.uniform(-8, 8))) for _ in range(3)]
    atoms = []
    for el, *r in molecule:
        r = [D(v) for v in r]
        xyz = [sum(rot[i][k] * r[k] for k in range(3)) + shift[i] for i in range(3)]
        if decimals:
            xyz = [v.quantize(D(1).scaleb(-decimals), rounding=ROUND_HALF_UP) for v in xyz]
        atoms.append(Atom(el, *xyz, CHARGE[el]))
    rng.shuffle(atoms)
    getcontext().prec = 6
    return atoms


def same_atom_set(ref, out, tol):
    unused = list(out)
    for r in ref:
        dist, o = min((max(abs(r.x - o.x), abs(r.y - o.y), abs(r.z - o.z)), i)
                      for i, o in enumerate(unused) if o.element == r.element)
        if dist > tol:
            return False
        unused.pop(o)
    return not unused


def test_sno_rotation_invariance():
    orig_prec = getcontext().prec
    rng = random.Random(0)
    try:
        for name, molecule in MOLECULES.items():
            for decimals in (None, 7):
                outs = [orient_system(moved(molecule, rng, decimals)) for _ in range(10)]
                if decimals is None:
                    strings = {tuple((a.element, str(a.x), str(a.y), str(a.z)) for a in o) for o in outs}
                    assert len(strings) == 1, f"{name}: {len(strings)} distinct outputs"
                else:
                    for o in outs[1:]:
                        assert same_atom_set(outs[0], o, D("5e-6")), f"{name} (rounded input) not invariant"
    finally:
        getcontext().prec = orig_prec
