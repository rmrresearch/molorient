"""
Decimal-only translation of NWChem's automatic molecular orientation
(geom_auto_sym -> HND_AUTSYM -> HND_MOLOPS, mass-weighted). Citations are
geom_input.F / geom.F line numbers from nwchemgit/nwchem @ 76d77c1.

Deviations from a literal translation:
  - Axis order is found by rotate-and-compare instead of NWChem's complex
    power sums (no Decimal complex type).
  - The post-hoc fixup rotations use only the booleans needed to recognize
    the eight affected families, not HND_MOLSMB's full naming tree.
  - I/Ih get no special handling (NWChem's autosym errors on them).
"""

from decimal import Decimal, getcontext, ROUND_HALF_UP
import periodictable as pt

from molorient.classes.vector import Vector
from molorient.classes.square_matrix import SquareMatrix
from molorient.classes.atom import Atom
from molorient.utils.diagonalization import eigval_solver, eigvec_solver
from molorient.utils.translation import translation_vector, translate_to_origin
from molorient.utils.axis_standardization import orient_atom
from molorient.utils.trig_helpers import cos_series, sin_series, pi_as_decimal


# NWChem thresholds are in bohr / amu*bohr^2; rescaled to Angstrom with
# NWChem's angstrom_to_au (geom.F:206).
_ANG2AU = Decimal('1.88972598858')
THREQUIV = Decimal('1e-2') / _ANG2AU               # atom matching (159)
THREQUIV_MOMENT = Decimal('1e-2') / _ANG2AU ** 2   # symmetric-top unique axis (4962)
TENM04 = Decimal('1e-4') / _ANG2AU ** 2            # zero/degenerate moments (4553)
CONV = Decimal('1e-4') / _ANG2AU ** 2              # canon axes / handedness (6194)
TENM05_DIST = Decimal('1e-5') / _ANG2AU            # atom at origin (4976)
TENM05 = Decimal('1e-5')                           # off-axis unit component (4991)


# Most-abundant-isotope masses, geom.F:4694-4769 (not periodictable's
# standard atomic weights).
_DEF_MASSES = {
    1: '1.007825', 2: '4.0026', 3: '7.016', 4: '9.01218', 5: '11.00931',
    6: '12.0', 7: '14.00307', 8: '15.99491', 9: '18.9984', 10: '19.99244',
    11: '22.9898', 12: '23.98504', 13: '26.98154', 14: '27.97693', 15: '30.97376',
    16: '31.97207', 17: '34.96885', 18: '39.9624', 19: '38.96371', 20: '39.96259',
    21: '44.95592', 22: '45.948', 23: '50.9440', 24: '51.9405', 25: '54.9381',
    26: '55.9349', 27: '58.9332', 28: '57.9353', 29: '62.9298', 30: '63.9291',
    31: '68.9257', 32: '73.9219', 33: '74.9216', 34: '78.9183', 35: '79.9165',
    36: '83.912', 37: '84.9117', 38: '87.9056', 39: '88.9054', 40: '89.9043',
    41: '92.9060', 42: '97.9055', 43: '97.9072', 44: '101.9037', 45: '102.9048',
    46: '105.9032', 47: '106.90509', 48: '113.9036', 49: '114.9041', 50: '117.9018',
    51: '120.9038', 52: '129.9067', 53: '126.9004', 54: '131.9042', 55: '132.9051',
    56: '137.9050', 57: '138.9061', 58: '139.9053', 59: '140.9074', 60: '143.9099',
    61: '144.9128', 62: '151.9195', 63: '152.9209', 64: '157.9241', 65: '159.9250',
    66: '163.9288', 67: '164.9303', 68: '165.9304', 69: '168.9344', 70: '173.9390',
    71: '174.9409', 72: '179.9468', 73: '180.948', 74: '183.9510', 75: '186.9560',
    76: '189.9586', 77: '192.9633', 78: '194.9648', 79: '196.9666', 80: '201.9706',
    81: '204.9745', 82: '207.9766', 83: '208.9804', 84: '209.9829', 85: '210.9875',
    86: '222.0175', 87: '223.0198', 88: '226.0254', 89: '227.0278', 90: '232.0382',
    91: '231.0359', 92: '238.0508', 93: '237.0482', 94: '244.0642', 95: '243.0614',
    96: '247.0704', 97: '247.0703', 98: '251.0796', 99: '252.0829', 100: '257.0950',
    101: '258.0986', 102: '259.1009', 103: '262.1100', 104: '261.1087', 105: '262.1138',
    106: '266.1219', 107: '262.1229', 108: '267.1318', 109: '268.1388',
    110: '281.0000', 111: '280.0000', 112: '285.0000', 113: '284.0000', 114: '289.0000',
    115: '287.0000', 116: '292.0000', 117: '293.0000', 118: '294.0000',
    119: '0', 120: '0',
}


def element_mass(atom):
    """NWChem default mass (geom.F:4694-4769); ghost centers get 0."""
    atn = pt.elements.symbol(atom.element).number
    return Decimal(_DEF_MASSES[atn]) if atn else Decimal('0')


def _vec(elements):
    v = Vector(3)
    v.elements = list(elements)
    return v


def _pos_vector(atom):
    return _vec((atom.x, atom.y, atom.z))


def _rotate_vector(pos_vec, rot_mat):
    return rot_mat.transpose().multiply(pos_vec)


def _rot_mat_from_axes(axes):
    """axes[0], axes[1], axes[2] become the X, Y, Z columns."""
    rot_mat = SquareMatrix(3)
    rot_mat.elements = [[axes[col].elements[row] for col in range(3)] for row in range(3)]
    return rot_mat


def _negate(pos, k):
    return _vec(-c if i == k else c for i, c in enumerate(pos.elements))


def _maps_onto(images, positions, elements, tol):
    """True if every atom's image lands within tol of a same-element atom."""
    return all(
        any(elements[j] == elements[i]
            and all(abs(a - b) < tol for a, b in zip(image.elements, positions[j].elements))
            for j in range(len(positions)))
        for i, image in enumerate(images))


def _mass_centered(atoms):
    """Shift to the center of mass for analysis (geom_input.F:4375-4378)."""
    getcontext().prec += 5
    masses = [element_mass(a) for a in atoms]
    total = sum(masses)
    cx, cy, cz = (sum(m * getattr(a, c) for m, a in zip(masses, atoms)) / total for c in 'xyz')
    shifted = [Atom(a.element, a.x - cx, a.y - cy, a.z - cz, a.charge) for a in atoms]
    getcontext().prec -= 5
    return shifted


def build_tensor(atoms):
    """Mass-weighted inertia tensor (geom.F:7527-7591)."""
    masses = [element_mass(a) for a in atoms]
    I_xx = sum(m * (a.y**2 + a.z**2) for m, a in zip(masses, atoms))
    I_yy = sum(m * (a.x**2 + a.z**2) for m, a in zip(masses, atoms))
    I_zz = sum(m * (a.x**2 + a.y**2) for m, a in zip(masses, atoms))
    I_xy = -sum(m * a.x * a.y for m, a in zip(masses, atoms))
    I_xz = -sum(m * a.x * a.z for m, a in zip(masses, atoms))
    I_yz = -sum(m * a.y * a.z for m, a in zip(masses, atoms))

    tensor = SquareMatrix(3)
    tensor.elements = [[I_xx, I_xy, I_xz], [I_xy, I_yy, I_yz], [I_xz, I_yz, I_zz]]
    return tensor


def _swap(eigvals, eigvecs, i, j):
    eigvals[i], eigvals[j] = eigvals[j], eigvals[i]
    eigvecs[i], eigvecs[j] = eigvecs[j], eigvecs[i]


def canon_axes(eigvals, eigvecs):
    """geom_canon_axes (geom_input.F:5796-5843): fix signs, order degenerate axes by nodes."""
    getcontext().prec += 5

    nodes = []
    for i in range(3):
        a = eigvecs[i].elements
        nodes.append((a[0] * a[1] < 0) + (a[1] * a[2] < 0))
        if 3 * a[0] + 2 * a[1] + a[2] < 0:
            eigvecs[i] = eigvecs[i].scale(-1)

    for i in range(2):
        for j in range(i + 1, 3):
            if abs(eigvals[j] - eigvals[i]) < CONV and nodes[j] < nodes[i]:
                _swap(eigvals, eigvecs, i, j)

    getcontext().prec -= 5
    return eigvals, eigvecs


def _diagonalize(tensor):
    """
    HND_MOLAXS (geom_input.F:6174-6254): eigenvalues DESCENDING (6207-6220),
    canon_axes, then force a right-handed frame (6223-6252).
    """
    eigvals = sorted(eigval_solver(tensor), reverse=True)
    eigvecs = list(eigvec_solver(*eigvals, tensor))
    eigvals, eigvecs = canon_axes(eigvals, eigvecs)

    getcontext().prec += 5
    if eigvecs[0].dot(eigvecs[1].cross(eigvecs[2])) <= 0:
        if abs(eigvals[0] - eigvals[1]) <= CONV:
            _swap(eigvals, eigvecs, 0, 1)
        elif abs(eigvals[1] - eigvals[2]) <= CONV:
            _swap(eigvals, eigvecs, 1, 2)
        else:
            eigvecs[2] = eigvecs[2].scale(-1)
    getcontext().prec -= 5
    return eigvals, eigvecs


def mass_inertia_tensor(atoms):
    getcontext().prec += 10
    eigvals, eigvecs = _diagonalize(build_tensor(atoms))
    getcontext().prec -= 10
    return eigvals, eigvecs


def classify_top(eigvals):
    """geom_input.F:4595-4623. Linear molecules also pass as 'symmetric' downstream."""
    getcontext().prec += 4
    nzer = sum(abs(e) < TENM04 for e in eigvals)
    deg01 = abs(eigvals[1] - eigvals[0]) < TENM04
    deg12 = abs(eigvals[1] - eigvals[2]) < TENM04
    getcontext().prec -= 4

    if nzer > 1:
        return 'atom'
    if nzer == 1:
        return 'linear'
    if deg01 and deg12:
        return 'spherical'
    return 'symmetric' if deg01 or deg12 else 'asymmetric'


def orient_symmetric_or_linear(eigvals, eigvecs, atoms):
    """
    geom_input.F:4954-5026. Unique axis K from the degenerate pair; Y from the
    first off-axis atom; the remaining axis is Y x K (negated when K is X).
    Linear molecules find no off-axis atom and keep their eigenvectors.
    """
    getcontext().prec += 5

    kaxis = None
    if abs(eigvals[1] - eigvals[0]) < THREQUIV_MOMENT:
        kaxis = 2
    if abs(eigvals[1] - eigvals[2]) < THREQUIV_MOMENT:
        kaxis = 0
    if kaxis is None:
        getcontext().prec -= 5
        raise ValueError("orient_symmetric_or_linear: no degenerate pair found")

    kvec = eigvecs[kaxis]
    other = 2 - kaxis
    axes = list(eigvecs)
    for atom in atoms:
        pos = _pos_vector(atom)
        norm = pos.dot(pos).sqrt()
        if norm <= TENM05_DIST:
            continue
        unit = pos.scale(1 / norm)
        perp = unit.add(kvec.scale(-unit.dot(kvec)))
        perp_norm = perp.dot(perp).sqrt()
        if perp_norm > TENM05:
            axes[1] = perp.scale(1 / perp_norm)
            cross = axes[1].cross(kvec)
            axes[other] = cross if kaxis == 2 else cross.scale(-1)
            break

    rot_mat = _rot_mat_from_axes(axes)
    getcontext().prec -= 5
    return rot_mat


def orient_spherical(eigvals, eigvecs, atoms):
    """
    Cubic groups, geom_input.F:4625-4951. Push each equivalent atom pair 1%
    outward, rediagonalize, and keep candidate axes that are real C2/C4/S4
    axes of the original molecule; two confirmed axes fix the frame.
    Returns None when no valid frame is found (4938-4943).
    """
    getcontext().prec += 5
    tol = THREQUIV

    n = len(atoms)
    positions = [_pos_vector(a) for a in atoms]
    elements = [a.element for a in atoms]
    distances = [p.dot(p).sqrt() for p in positions]
    equivalence = list(range(n))
    for i in range(n):
        if equivalence[i] != i:
            continue
        if distances[i] <= tol:
            equivalence[i] = -1
            continue
        for j in range(n):
            if j != i and elements[j] == elements[i] and abs(distances[j] - distances[i]) < tol:
                equivalence[j] = i

    syminv = _maps_onto([p.scale(-1) for p in positions], positions, elements, tol)

    def frame(ax, a, j, k):
        e = [None] * 3
        e[ax], e[(ax + 1) % 3], e[(ax + 2) % 3] = a, j, k
        return _vec(e)

    def axis_tests(axes):
        rot_mat = _rot_mat_from_axes(axes)
        rotated = [_rotate_vector(p, rot_mat) for p in positions]
        xs = [p.elements for p in rotated]
        c2, c4, s4 = [], [], []
        for ax in range(3):
            J, K = (ax + 1) % 3, (ax + 2) % 3
            c2.append(_maps_onto([frame(ax, x[ax], -x[J], -x[K]) for x in xs], rotated, elements, tol))
            c4.append(_maps_onto([frame(ax, x[ax], -x[K], x[J]) for x in xs], rotated, elements, tol))
            s4.append(_maps_onto([frame(ax, -x[ax], -x[K], x[J]) for x in xs], rotated, elements, tol))
        return c2, c4, s4

    pairs = ((i, j) for i in range(n) if equivalence[i] == i
             for j in range(n) if j != i and equivalence[j] == i)
    axm = []
    for i, j in pairs:
        s = Decimal('1.01')
        distorted = [Atom(a.element, a.x * s, a.y * s, a.z * s, a.charge) if k in (i, j) else a
                     for k, a in enumerate(atoms)]
        _, axes = _diagonalize(build_tensor(distorted))
        c2, c4, s4 = axis_tests(axes)

        if syminv:  # T_h / O_h
            keep = [c2[k] and (c4[k] or not any(c4)) for k in range(3)]
        else:       # T / T_d / O
            grpo, grptd = any(c4), any(s4)
            keep = [c2[k] and ((grpo and c4[k]) or (grptd and s4[k]) or not (grpo or grptd))
                    for k in range(3)]
        axm += [axes[k] for k in range(3) if keep[k]]
        if len(axm) >= 2:
            break

    getcontext().prec -= 5

    if len(axm) < 2:
        return None

    axis0, axis1 = axm[:2]
    axis2 = axis0.cross(axis1)
    det = axis0.dot(axis1.cross(axis2))
    if det <= Decimal('0.99') or det > Decimal('1.1'):
        return None

    return _rot_mat_from_axes([axis0, axis1, axis2])


def _quantized_atoms(atoms, coords):
    """Atoms at coords, rounded to the caller's precision (call at prec + 10)."""
    tol = Decimal(1).scaleb(-(getcontext().prec - 10))
    return [Atom(a.element, *(c.quantize(tol, rounding=ROUND_HALF_UP) for c in xyz), a.charge)
            for a, xyz in zip(atoms, coords)]


def rotate_atoms(atoms, rot_mat):
    """Rotate atoms into the frame whose axes are rot_mat's columns."""
    getcontext().prec += 10
    rotated = _quantized_atoms(
        atoms, [_rotate_vector(_pos_vector(a), rot_mat).elements for a in atoms])
    getcontext().prec -= 10
    return rotated


def _rotate_in_plane(positions, axis_index, theta):
    """Rotate by theta about axis_index (rot_theta_z, geom_input.F:7193-7212, any axis)."""
    i, j = [k for k in range(3) if k != axis_index]
    cos_t, sin_t = cos_series(theta), sin_series(theta)
    rotated = []
    for p in positions:
        e = list(p.elements)
        e[i] = cos_t * p.elements[i] - sin_t * p.elements[j]
        e[j] = sin_t * p.elements[i] + cos_t * p.elements[j]
        rotated.append(_vec(e))
    return rotated


def detect_axis_order(positions, elements, axis_index, tol, max_order=24):
    """
    Highest n <= max_order whose 2*pi/n rotation (or rotation-reflection)
    about axis_index maps the molecule onto itself; returns (order, proper).
    Replaces the complex power sums of geom_input.F:5061-5226, so on-axis
    atoms need no linear-molecule override (5127-5144).
    """
    best_order, best_proper = 1, True
    for order in range(2, max_order + 1):
        rotated = _rotate_in_plane(positions, axis_index, 2 * pi_as_decimal() / order)
        proper = _maps_onto(rotated, positions, elements, tol)
        improper = _maps_onto([_negate(p, axis_index) for p in rotated], positions, elements, tol)
        if proper or improper:
            best_order, best_proper = order, proper
    return best_order, best_proper


def align_highest_order_to_z(orders):
    """geom_input.F:5300-5336: permutation putting the highest-order axis on Z (ties -> X)."""
    if orders[2] > orders[0] and orders[2] > orders[1]:
        kaxis = 2
    elif orders[1] > orders[0] and orders[1] > orders[2]:
        kaxis = 1
    else:
        kaxis = 0
    return [(kaxis + 1 + k) % 3 for k in range(3)]


def _permute_positions(positions, newaxs):
    return [_vec(p.elements[k] for k in newaxs) for p in positions]


def find_perpendicular_c2(positions, elements, orders, propers, tol):
    """
    geom_input.F:5436-5489: X or Y as an even-order proper axis, else scan
    about Z in 0.5-degree steps for a perpendicular C2 (check_c2_perp, 7213).
    """
    symc2x = orders[0] % 2 == 0 and propers[0]
    symc2y = orders[1] % 2 == 0 and propers[1]
    if symc2x or symc2y or orders[2] <= 1:
        return symc2x, symc2y, positions

    step = Decimal('0.5') * (2 * pi_as_decimal() / 360)
    n_steps = int((2 * pi_as_decimal() / orders[2]) / step)
    for istep in range(n_steps + 1):
        trial = _rotate_in_plane(positions, 2, istep * step)
        for axis in (0, 1):
            if _maps_onto(_rotate_in_plane(trial, axis, pi_as_decimal()), trial, elements, tol):
                return axis == 0, axis == 1, trial
    return False, False, positions


def _cs_correction(positions, mirryz, mirrzx, mirrxy, order_z):
    """geom_input.F:5622-5648: with no rotation axis, move a YZ or ZX mirror onto XY."""
    if order_z != 1 or not (mirryz or mirrzx):
        return positions, mirryz, mirrzx, mirrxy
    if mirryz:
        return _permute_positions(positions, [2, 1, 0]), False, mirrzx, True
    return _permute_positions(positions, [0, 2, 1]), mirryz, False, True


def apply_posthoc_fixup(positions, order_z, proper_z, mirryz, mirrzx, mirrxy, symc2x, symc2y):
    """
    geom_input.F:5667-5715: extra Z rotation for C3v/C6v/D3/D6, C5v, D4d,
    D5h/D5 and D5d. Dnd halves the order (HND_MOLSMB, 6005-6007). Least
    verified stage of this module.
    """
    is_dn = symc2x or symc2y
    is_dnd = is_dn and (mirryz or mirrzx or not proper_z) and not mirrxy
    n = order_z // 2 if is_dnd else order_z
    pi = pi_as_decimal()

    theta = None
    if n in (3, 6) and not mirrxy and not is_dnd and (mirryz or mirrzx or is_dn):
        theta = pi / 4 if symc2x or (not symc2y and mirrzx) else -pi / 4  # C3v/C6v/D3/D6
    elif n == 5 and not is_dn and mirryz:
        theta = 2 * pi / 4  # C5v
    elif n == 4 and is_dnd:
        theta = Decimal('22.5') / 180 * pi  # D4d
    elif n == 5 and is_dn and not is_dnd and symc2x:
        theta = 2 * pi / 4  # D5h / D5
    elif n == 5 and is_dnd and symc2y:
        theta = -(2 * pi / 4)  # D5d

    return positions if theta is None else _rotate_in_plane(positions, 2, theta)


def orient_mol(atoms):
    """
    Mass-centered analysis -> principal axes -> top-specific frame -> highest
    order axis on Z -> perpendicular C2 and mirrors -> Cs and family fixups ->
    re-center on nuclear charge (translation.py).
    """
    analysis_atoms = _mass_centered(atoms)
    eigvals, eigvecs = mass_inertia_tensor(analysis_atoms)
    top = classify_top(eigvals)

    if top == 'atom':
        return rotate_atoms(translate_to_origin(atoms, translation_vector(atoms)), orient_atom())

    rot_mat = None
    if top in ('linear', 'symmetric'):
        rot_mat = orient_symmetric_or_linear(eigvals, eigvecs, analysis_atoms)
    elif top == 'spherical':
        rot_mat = orient_spherical(eigvals, eigvecs, analysis_atoms)
    candidate_atoms = rotate_atoms(analysis_atoms, rot_mat or _rot_mat_from_axes(eigvecs))

    getcontext().prec += 5
    tol = THREQUIV

    positions = [_pos_vector(a) for a in candidate_atoms]
    elements = [a.element for a in candidate_atoms]

    orders, propers = zip(*(detect_axis_order(positions, elements, k, tol) for k in range(3)))

    newaxs = align_highest_order_to_z(orders)
    positions = _permute_positions(positions, newaxs)
    orders = [orders[k] for k in newaxs]
    propers = [propers[k] for k in newaxs]

    symc2x, symc2y, positions = find_perpendicular_c2(positions, elements, orders, propers, tol)
    mirryz, mirrzx, mirrxy = (_maps_onto([_negate(p, k) for p in positions], positions, elements, tol)
                              for k in range(3))
    positions, mirryz, mirrzx, mirrxy = _cs_correction(positions, mirryz, mirrzx, mirrxy, orders[2])
    positions = apply_posthoc_fixup(positions, orders[2], propers[2],
                                    mirryz, mirrzx, mirrxy, symc2x, symc2y)

    oriented = [Atom(a.element, *p.elements, a.charge) for a, p in zip(candidate_atoms, positions)]
    oriented = translate_to_origin(oriented, translation_vector(oriented))

    getcontext().prec += 5
    final_atoms = _quantized_atoms(oriented, [(a.x, a.y, a.z) for a in oriented])
    getcontext().prec -= 10
    return final_atoms