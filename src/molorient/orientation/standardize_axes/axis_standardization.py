from decimal import Decimal, getcontext, ROUND_HALF_UP
from functools import cmp_to_key
from itertools import combinations, permutations
from molorient.classes.atom import Atom
from molorient.classes.square_matrix import SquareMatrix
from molorient.classes.vector import Vector
from molorient.helpers.diagonalization import (
    eigval_solver, eigvec_solver, pos_vector, build_inertia_tensor, rot_mat_from_axes,
)
from molorient.helpers.trig_helpers import sin_series, cos_series, pi_as_decimal
from molorient.helpers.precision import prec_tol, moment_tol, geom_tol


def classify_moments(raw_moments, p):
    """
    Classifies the top from the unrounded principal moments (ascending) with a
    tolerance tied to the user's precision p. Returns (label, unique), where
    label is 'atom', 'linear', 'symmetric', 'spherical' or 'asymmetric' and
    unique is the index (0 or 2) of the non-degenerate moment of a symmetric
    top (0 for a linear molecule, None otherwise).
    """
    a, b, c = raw_moments
    if c == 0:
        return 'atom', None
    tau = moment_tol(p) * c
    if a <= tau and c - b <= tau:
        return 'linear', 0
    if c - a <= tau:
        return 'spherical', None
    if b - a <= tau:
        return 'symmetric', 2
    if c - b <= tau:
        return 'symmetric', 0
    return 'asymmetric', None


def snap_moments(raw_moments, label):
    """Raw moments with every degenerate group replaced by its exact mean."""
    a, b, c = raw_moments
    if label == 'atom':
        return [Decimal(0)] * 3
    if label == 'linear':
        mean = (b + c) / 2
        return [Decimal(0), mean, mean]
    if label == 'spherical':
        mean = (a + b + c) / 3
        return [mean] * 3
    if label == 'symmetric':
        if b - a <= c - b:
            mean = (a + b) / 2
            return [mean, mean, c]
        mean = (b + c) / 2
        return [a, mean, mean]
    return [a, b, c]


def unique_axis_frame(tensor, raw_moments, unique):
    """
    Right-handed orthonormal axes [e0, e1, e2] with the unique axis of a linear
    or symmetric top in column `unique` (0 or 2). The unique axis comes from
    the projector (A - l1)(A - l2) onto it, built from the two degenerate
    moments, so it stays accurate when they are only equal within noise.
    """
    others = [m for i, m in enumerate(raw_moments) if i != unique]
    m_1 = [[tensor.elements[i][j] - (others[0] if i == j else 0) for j in range(3)] for i in range(3)]
    m_2 = [[tensor.elements[i][j] - (others[1] if i == j else 0) for j in range(3)] for i in range(3)]
    columns = []
    for j in range(3):
        col = Vector(3)
        for i in range(3):
            col.elements[i] = sum(m_1[i][k] * m_2[k][j] for k in range(3))
        columns.append(col)
    Z = max(columns, key=lambda c: c.dot(c)).normalized()

    basis = []
    for i in range(3):
        e = Vector(3)
        e.assign(i, 1)
        basis.append(e)
    e = min(basis, key=lambda v: abs(Z.dot(v)))
    X = e.subtract(Z.scale(Z.dot(e))).normalized()
    Y = Z.cross(X)

    return [X, Y, Z] if unique == 2 else [Z, X, Y]


def inertia_tensor(atoms, p=None):
    """
    This function calculates the inertia tensor of the system using charge instead of mass. It returns
    the principal moments of inertia (eigenvalues) of the tensor in order of increasing magnitude.
    Moments equal to within a tolerance tied to the user precision p (default: the current Decimal
    precision) are returned exactly equal, so downstream classification can compare them directly.
    """

    if p is None:
        p = getcontext().prec
    atoms = sorted(atoms, key=lambda a: (a.element, a.x, a.y, a.z))
    getcontext().prec += 10

    tensor = build_inertia_tensor(atoms, [atom.charge for atom in atoms])

    raw_eigvals = list(sorted(eigval_solver(tensor)))
    label, unique = classify_moments(raw_eigvals, p)
    eigvals = snap_moments(raw_eigvals, label)

    if label in ('linear', 'symmetric'):
        eigvecs = unique_axis_frame(tensor, raw_eigvals, unique)
    else:
        eigvecs = list(eigvec_solver(*raw_eigvals, tensor))

    getcontext().prec -= 10

    #Round values
    #(1) If number of decimal places > number of sig figs, round to 10**(-precision)
    #(2) If number of sig figs > number of decimal places, round to nearest precision
    for i, e in enumerate(eigvals):
        t = e.as_tuple()
        sig_figs = len(t.digits)
        dec_places = max(0, -t.exponent)
        if sig_figs < dec_places:
            eigvals[i] = e.quantize(prec_tol())
        else:
            rounded = round(e, getcontext().prec - e.adjusted() - 1)
            eigvals[i] = Decimal(str(rounded))

    #Assign eigenvalues to their correspoinding eigenvector for asymmetric top.
    if label == 'asymmetric':
        remaining = list(eigvecs)
        ordered = []

        for e in raw_eigvals:
            def residual(v):
                Av = tensor.multiply(v)
                return sum(abs(Av.elements[k] - e * v.elements[k]) for k in range(3))

            best_index = min(
                range(len(remaining)),
                key=lambda i: residual(remaining[i])
            )

            best = remaining.pop(best_index)
            ordered.append(best)

        eigvecs = list(ordered)

    #Assign right-handed eigenbasis coordinate system
    det = (eigvecs[0].cross(eigvecs[1])).dot(eigvecs[2])
    if det < 0:
        eigvecs[2] = eigvecs[2].scale(-1)

    return eigvals, eigvecs


def orient_atom():
    """
    If all moments == 0, the system is a single atom. No rotation applied.
    """
    return SquareMatrix.identity(3)


def _tolerant_compare(a, b, tol):
    d = a - b
    if abs(d) <= tol:
        return 0
    return -1 if d < 0 else 1


def _frame_key(frame, atoms, tol):
    """The atoms rotated into frame, sorted by (charge, x, y, z) with ties within tol."""
    rows = []
    rot = frame.transpose()
    for atom in atoms:
        new = rot.multiply(pos_vector(atom)).elements
        rows.append((atom.charge, new[0], new[1], new[2]))

    def compare_rows(r, s):
        if r[0] != s[0]:
            return -1 if r[0] < s[0] else 1
        for u, v in zip(r[1:], s[1:]):
            c = _tolerant_compare(u, v, tol)
            if c:
                return c
        return 0

    rows.sort(key=cmp_to_key(compare_rows))
    return rows, compare_rows


def pick_canonical_frame(frames, atoms, tol):
    """
    The frame, out of equally valid candidates, that gives the lexicographically
    smallest atom list (atoms sorted by charge, x, y, z; coordinates within tol
    count as equal). Symmetry-equivalent candidates give the same geometry, so
    the choice does not depend on how the molecule was originally oriented.
    Every candidate must be a proper rotation.
    """
    best_frame = None
    best_rows = None
    for frame in frames:
        cols = [Vector(3) for _ in range(3)]
        for i in range(3):
            for k in range(3):
                cols[i].elements[k] = frame.elements[k][i]
        if cols[0].cross(cols[1]).dot(cols[2]) <= 0:
            raise ValueError("pick_canonical_frame: candidate frame is not a proper rotation.")
        rows, compare_rows = _frame_key(frame, atoms, tol)
        if best_rows is None:
            best_frame, best_rows = frame, rows
            continue
        for r, s in zip(rows, best_rows):
            c = compare_rows(r, s)
            if c:
                if c < 0:
                    best_frame, best_rows = frame, rows
                break

    return best_frame


def orient_linear(eigvecs, atoms, p=None):
    """
    If moment_a == 0 and moment_b == moment_c, the system is linear. The vector through
    the line of atoms (eigvecs[0]) acts as the Z axis; its direction is whichever of the two
    choices gives the smaller sorted geometry.
    """
    if p is None:
        p = getcontext().prec
    tol = geom_tol(atoms, p)
    Z, X, Y = eigvecs[0], eigvecs[1], eigvecs[2]
    frames = [rot_mat_from_axes([X, Y, Z]), rot_mat_from_axes([X, Y.scale(-1), Z.scale(-1)])]

    return pick_canonical_frame(frames, atoms, tol)


def symmetric_frames(Z, atoms, tol):
    """
    Every candidate frame for a symmetric top with axis Z. The Y axis points at an
    off-axis atom chosen by: (1) nearest XY plane, (2) positive Z projection,
    (3) nearest Z axis, (4) lowest atomic number. Z's sign is arbitrary, so both
    senses are tried and a sense is kept only if criterion (2) can be met (both
    are kept if neither can). X = Y cross Z. Returns one frame [X, Y, Z] per
    candidate atom and kept sense; the caller picks among them.
    """
    info = []
    for atom in atoms:
        pos = pos_vector(atom)
        z_proj = Z.dot(pos)
        perp = pos.subtract(Z.scale(z_proj))
        dist_to_z = perp.norm()
        if dist_to_z > tol:
            info.append((atom, z_proj, dist_to_z, perp))
    if not info:
        raise ValueError("symmetric_frames: every atom lies on the symmetry axis.")

    def candidates(sense):
        #(1) Nearest to XY plane.
        min_d_xy = min(abs(t[1]) for t in info)
        found = [t for t in info if abs(t[1]) - min_d_xy <= tol]

        #(2) Positive Z projection.
        positive = [t for t in found if sense * t[1] > tol]
        if positive:
            found = positive

        #(3) Nearest to Z axis
        min_dz = min(t[2] for t in found)
        found = [t for t in found if t[2] - min_dz <= tol]

        #(4) Lowest atomic number
        min_charge = min(t[0].charge for t in found)
        found = [t for t in found if t[0].charge == min_charge]

        return found, bool(positive)

    results = {sense: candidates(sense) for sense in (1, -1)}
    senses = [s for s in (1, -1) if results[s][1]] or [1, -1]

    frames = []
    for sense in senses:
        Z_s = Z.scale(sense)
        for atom, z_proj, dist_to_z, perp in results[sense][0]:
            Y = perp.normalized()
            X = Y.cross(Z_s).normalized()
            frames.append(rot_mat_from_axes([X, Y, Z_s]))

    return frames


def orient_symm(moment_a, moment_b, eigvecs, atoms, p=None):
    """
    If two moments are equal and one is unequal, the unequal moment acts as the Z axis. The Y axis is determined by the following process:
    Like atoms lying perpendicular to the Z axis are candidates for the Y axis. To remove ambiguity of the Y axis, the following criteria
    are imposed upon the candidates until a winner is chosen:
    (1) Nearest to XY plane.
    (2) Positive Z axis projection.
    (3) Nearest to Z axis
    (4) Lowest atomic number.
    The winning atom serves as the direction of the Y axis (all remaining ties, and the unphysical sign of Z, are
    settled by pick_canonical_frame). The X axis is simply a cross product of Y and Z.
    """
    if p is None:
        p = getcontext().prec
    getcontext().prec += 5
    tol = geom_tol(atoms, p)
    z_col = 2 if moment_a == moment_b else 0
    Z = Vector(3)
    for i in range(3):
        Z.elements[i] = eigvecs[z_col].elements[i]

    rot_mat = pick_canonical_frame(symmetric_frames(Z, atoms, tol), atoms, tol)

    getcontext().prec -= 5

    return rot_mat


def orient_spherical(atoms, group, axes, p=None):
    """
    If all moments are equal, the system is of high symmetry. The principal axes of rotation are used.
    For tetrahedral symmetry, the 3 C2 axes are used as the Cartesian axes. For octahedral symmetry, the 3 C4
    axes are used. For icosahedral symmetry, a C5 axis is chosen as the Z axis and the system is treated
    as a symmetric top. Among all equivalent choices, pick_canonical_frame selects one, so the result is
    independent of the input orientation (and always a proper rotation).
    """
    if p is None:
        p = getcontext().prec
    getcontext().prec += 2
    try:
        tol = geom_tol(atoms, p)

        #Td/Oh point groups: the three C2 (Td) or C4 (Oh) axes are used as the Cartesian axes
        if group in ('Td', 'Oh'):
            frames = []
            for i, j in permutations(range(3), 2):
                for sign_x in (1, -1):
                    for sign_y in (1, -1):
                        X = axes[i].scale(sign_x)
                        Y = axes[j].scale(sign_y)
                        Y = Y.subtract(X.scale(X.dot(Y))).normalized()
                        frames.append(rot_mat_from_axes([X, Y, X.cross(Y)]))

        #Ih point group: A C5 axis is chosen and used as Z. The system
        #is treated as a symmetric top.
        elif group == 'Ih':
            frames = []
            for Z in axes:
                frames.extend(symmetric_frames(Z, atoms, tol))

        else:
            raise ValueError(f"orient_spherical: unknown point group {group!r}.")

        rot_mat = pick_canonical_frame(frames, atoms, tol)
    finally:
        getcontext().prec -= 2

    return rot_mat


def orient_asymm(eigvecs, atoms, p=None):
    """
    If all moments are unequal, the principal axes are used as the Cartesian axes. Each axis
    has an arbitrary sign; of the four proper sign choices, pick_canonical_frame selects the
    one giving the smallest sorted geometry, comparing the unrounded coordinates within
    geom_tol (so noise cannot decide between coordinates that are equal by symmetry).
    """
    if p is None:
        p = getcontext().prec
    getcontext().prec += 2
    try:
        tol = geom_tol(atoms, p)
        frames = []
        for signs in ((1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)):
            frames.append(rot_mat_from_axes([v.scale(s) for v, s in zip(eigvecs, signs)]))
        return pick_canonical_frame(frames, atoms, tol)
    finally:
        getcontext().prec -= 2


def standardize_axes(moments, eigvecs, atoms, p=None):
    """
    This function uses the principal moments of nuclear inertia to classify the type
    of top of the system: asymmetric, symmetric, spherical, linear, or a single
    atom. The system is then rotated to align the Cartesian axes with the principal
    axes.
    """

    if p is None:
        p = getcontext().prec

    moment_a = moments[0]
    moment_b = moments[1]
    moment_c = moments[2]

    #Single atom
    if moment_a == moment_b == moment_c == 0:
        rot_mat = orient_atom()
    
    #Linear
    elif (moment_a == 0) and moment_b == moment_c:
        rot_mat = orient_linear(eigvecs, atoms, p)
    
    #Symmetric top
    elif (moment_a == moment_b != moment_c) or (moment_a != moment_b == moment_c):
        rot_mat = orient_symm(moment_a, moment_b, eigvecs, atoms, p)
    
    #Spherical top
    elif moment_a == moment_b == moment_c:
        group, axes = cn_axes_finder(atoms, p)
        rot_mat = orient_spherical(atoms, group, axes, p)
    
    #Asymmetric top
    elif moment_a != moment_b != moment_c:
        rot_mat = orient_asymm(eigvecs, atoms, p)

    else:
        raise ValueError(f"standardize_axes: moments {moments} did not match any top classification.")

    #Rotation
    standardized_atoms = []
    getcontext().prec += 10
    tol = prec_tol(10)

    for atom in atoms:
        pos_vec = pos_vector(atom)

        new_pos = (rot_mat.transpose()).multiply(pos_vec)

        standardized_atoms.append(Atom(atom.element, 
                                        new_pos.elements[0].quantize(tol, rounding = ROUND_HALF_UP),
                                        new_pos.elements[1].quantize(tol, rounding = ROUND_HALF_UP),
                                        new_pos.elements[2].quantize(tol, rounding = ROUND_HALF_UP),
                                        atom.charge))

    getcontext().prec -= 10
    
    return standardized_atoms


def find_cn_axes(candidates, atoms, order, tol):
    """
    Which of candidates are true order-fold rotation axes of atoms, tested
    via Rodrigues' rotation formula. An atom's rotated image must land
    within distance sqrt(tol) of an atom (tol is a squared distance).
    """
    theta = 2 * pi_as_decimal() / order
    cos_theta = cos_series(theta)
    sin_theta = sin_series(theta)
    found = []
    for k in candidates:
        is_axis = True
        for atom in atoms:
            v = pos_vector(atom)
            v_rot = v.scale(cos_theta).add((k.cross(v)).scale(sin_theta).add(k.scale((k.dot(v)) * (1 - cos_theta))))

            if not any(
                (v_rot.elements[0] - atom2.x)**2 +
                (v_rot.elements[1] - atom2.y)**2 +
                (v_rot.elements[2] - atom2.z)**2 < tol
                for atom2 in atoms
            ):
                is_axis = False
                break

        if is_axis:
            found.append(k)
    return found


def cn_axes_finder(atoms, p=None):
    """
    This function is designed for the spherical top case. Finds principal rotation axes of high symmetry point groups.
    Identifies the system as tetrahedral, octahedral, or icosahedral and returns the principal axes.
    Distances are compared within geom_tol(atoms, p); p defaults to the current Decimal precision.
    """
    if p is None:
        p = getcontext().prec

    getcontext().prec += 2
    try:
        dist_tol = geom_tol(atoms, p)
        angle_tol = moment_tol(p)

        #Principal axes can only be either going through atoms or bisecting two atoms.
        #Find axes going through atoms
        candidates = []

        for atom in atoms:
            v = pos_vector(atom)
            norm = v.norm()
            if norm > dist_tol:
                candidates.append(v.normalized())

        #Uniqueness check for axes through atoms
        unique = []
        for v in candidates:
            if not any(abs(abs(v.dot(u)) - 1) < angle_tol for u in unique):
                unique.append(v)
        candidates = unique

        #Finds bisecting axes
        bisect_candidates = list(candidates)
        for v_0, v_1 in combinations(bisect_candidates, r = 2):
            w = v_0.add(v_1)
            if w.norm() > angle_tol:
                candidates.append(w.normalized())

        #Another uniqueness check
        unique = []
        for v in candidates:
            if not any(abs(abs(v.dot(u)) - 1) < angle_tol for u in unique):
                unique.append(v)
        candidates = unique

        #Rotation around the candidate principal axes using Rodrigues' rotation formula.
        #Icosahedral symmetry utilizes the C5 axes, octahedral C4, and tetrahedral C2 for their rotations.
        sq_tol = dist_tol ** 2
        #Test for icosahedral:
        c5 = find_cn_axes(candidates, atoms, 5, sq_tol)
        if len(c5) > 0:
            return 'Ih', c5

        #Test for octahedral:
        c4 = find_cn_axes(candidates, atoms, 4, sq_tol)
        if len(c4) > 0:
            return 'Oh', c4

        #Test for tetrahedral:
        c2 = find_cn_axes(candidates, atoms, 2, sq_tol)
        if len(c2) > 0:
            return 'Td', c2
    finally:
        getcontext().prec -= 2

    raise ValueError("cn_axes_finder: no Td, Oh or Ih rotation axes found; the moments of inertia are degenerate but the geometry is not a spherical top.")
