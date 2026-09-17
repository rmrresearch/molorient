from decimal import Decimal, getcontext, ROUND_HALF_UP
from itertools import combinations
from molorient.classes.atom import Atom
from molorient.classes.square_matrix import SquareMatrix
from molorient.classes.vector import Vector
from molorient.utils.diagonalization import (
    eigval_solver, eigvec_solver, pos_vector, build_inertia_tensor, rot_mat_from_axes,
)
from molorient.utils.trig_helpers import sin_series, cos_series, arccos_series, arctan2, pi_as_decimal
from molorient.utils.precision import prec_tol


def find_eigenvector_for(tensor, eigvecs, e_unique, tol):
    """The eigvecs entry v satisfying tensor@v == e_unique*v (within tol), or None."""
    for v in eigvecs:
        Av = tensor.multiply(v)
        if all(abs(Av.elements[k] - e_unique * v.elements[k]) < tol for k in range(3)):
            return v
    return None


def inertia_tensor(atoms):
    """
    This function calculates the inertia tensor of the system using charge instead of mass. It returns
    the principal moments of inertia (eigenvalues) of the tensor in order of increasing magnitude.
    """

    atoms = sorted(atoms, key=lambda a: (a.element, a.x, a.y, a.z))
    getcontext().prec += 10

    tensor = build_inertia_tensor(atoms, [atom.charge for atom in atoms])

    moment_a, moment_b, moment_c = sorted(eigval_solver(tensor))
    v_0, v_1, v_2 = eigvec_solver(moment_a, moment_b, moment_c, tensor)
    eigvals = [moment_a, moment_b, moment_c]
    eigvecs = [v_0, v_1, v_2]
    raw_eigvals = [moment_a, moment_b, moment_c]

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
    
    #Assign the eigenvalues to its corresponding eigenvector by solving Av = λv for symmetric top
    if eigvals[0] == eigvals[1] != eigvals[2]:
        tol = prec_tol(2)
        unique_vec = find_eigenvector_for(tensor, eigvecs, eigvals[2], tol)
        remaining = [v for v in eigvecs if v is not unique_vec]
        eigvecs = [remaining[0], remaining[1], unique_vec]

    if (eigvals[0] != eigvals[1] == eigvals[2]) and eigvals[0] != 0:
        tol = prec_tol(2)
        unique_vec = find_eigenvector_for(tensor, eigvecs, eigvals[0], tol)
        remaining = [v for v in eigvecs if v is not unique_vec]
        eigvecs = [unique_vec, remaining[0], remaining[1]]

    #Assign eigenvalues to their correspoinding eigenvector for asymmetric top.
    if eigvals[0] != eigvals[1] != eigvals[2]:
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

        v_0, v_1, v_2 = ordered
        eigvecs = [v_0, v_1, v_2]
            
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

    
def orient_linear(atoms):
    """
    If moment_a == 0 and moment_b == moment_c, the system is linear. The vector through
    the line of atoms acts as the Z axis.
    """
    
    getcontext().prec += 2
    rot_mat = SquareMatrix(3)
    
    vec = Vector(3)
    vec.elements[0] = (atoms[-1].x - atoms[0].x)
    vec.elements[1] = (atoms[-1].y - atoms[0].y)
    vec.elements[2] = (atoms[-1].z - atoms[0].z)

    norm_vec = vec.normalized()

    theta = arccos_series(norm_vec.elements[2])
    phi = arctan2(norm_vec.elements[1], norm_vec.elements[0])
    rot_mat.assign(0, 0, (cos_series(theta + phi) + cos_series(theta - phi)) / 2)
    rot_mat.assign(0, 1, sin_series(phi))
    rot_mat.assign(0, 2, (sin_series(theta + phi) + sin_series(theta - phi)) / 2)
    rot_mat.assign(1, 0, (-sin_series(theta + phi) + sin_series(theta - phi)) / 2)
    rot_mat.assign(1, 1, cos_series(phi))
    rot_mat.assign(1, 2, (cos_series(theta + phi) - cos_series(theta - phi)) / 2)
    rot_mat.assign(2, 0, -sin_series(theta))
    rot_mat.assign(2, 1, Decimal('0'))
    rot_mat.assign(2, 2, cos_series(theta))
    
    getcontext().prec -= 2

    return rot_mat


def symmetric_axis_frame(Z, atoms, tol, quantize_group_key=False):
    """
    Given a chosen Z axis, groups atoms by (element, Z-projection) and picks a
    Y-axis candidate via: (1) nearest XY plane, (2) positive Z projection,
    (3) nearest Z axis, (4) lowest atomic number. X = Y_norm cross Z. Returns
    the rotation matrix with columns [X_norm, Y_norm, Z]. quantize_group_key
    matches orient_symm's original quantized grouping key; orient_spherical's
    Ih branch never quantized its key, so it passes False (kept as a
    parameter rather than unified, since the two are not the same value).
    """
    groups = {}
    for atom in atoms:
        pos = pos_vector(atom)
        Z_proj = Z.dot(pos)
        along = Z.scale(Z_proj)
        perp = pos.subtract(along)
        dist_to_z = (perp.dot(perp)).sqrt()
        if dist_to_z < tol:
            continue

        Z_key = Z_proj.quantize(tol) if quantize_group_key else Z_proj
        key = (atom.element, Z_key)
        if key not in groups:
            groups[key] = {
                'atoms': [],
                'Z_proj': Z_proj,
                'dist_to_xy': abs(Z_proj),
                'dist_to_Z': dist_to_z
            }
        groups[key]['atoms'].append(atom)

    candidates = list(groups.values())

    #(1) Nearest to XY plane.
    min_d_xy = min(g['dist_to_xy'] for g in candidates)
    candidates = [g for g in candidates if abs(g['dist_to_xy'] - min_d_xy) < tol]

    #(2) Positive Z projection.
    pos_z = [g for g in candidates if g['Z_proj'] > 0]
    if pos_z:
        candidates = pos_z

    #(3) Nearest to Z axis
    min_dz = min(g['dist_to_Z'] for g in candidates)
    candidates = [g for g in candidates if abs(g['dist_to_Z'] - min_dz) < tol]

    #(4) Lowest atomic number
    min_charge = min(atom.charge for g in candidates for atom in g['atoms'])
    candidates = [g for g in candidates if any(atom.charge == min_charge for atom in g['atoms'])]

    #Define Y axis
    key_atom = candidates[0]['atoms'][0]
    key_pos = pos_vector(key_atom)

    key_Z_proj = Z.dot(key_pos)
    Y = key_pos.add(Z.scale(-key_Z_proj))
    Y_norm = Y.normalized()
    X = Y_norm.cross(Z)
    X_norm = X.normalized()

    return rot_mat_from_axes([X_norm, Y_norm, Z])


def orient_symm(moment_a, moment_b, eigvecs, atoms):
    """
    If two moments are equal and one is unequal, the unequal moment acts as the Z axis. The Y axis is determined by the following process:
    Groups of like atoms lying perpendicular to the Z axis are candidates for the Y axis. To remove ambiguity of the Y axis, the following criteria
    are imposed upon the candidates until a winner is chosen:
    (1) Nearest to XY plane.
    (2) Positive Z axis projection.
    (3) Nearest to Z axis
    (4) Lowest atomic number.
    An arbitrary atom in the winning group is chosen and will serve as the direction of the Y axis. The X axis is simply a cross product of Y and Z.
    """
    getcontext().prec += 5
    tol = prec_tol(10)
    if moment_a == moment_b:
        z_col = 2
    else:
        z_col = 0
    Z = Vector(3)
    for i in range(3):
        Z.elements[i] = eigvecs[z_col].elements[i]

    rot_mat = symmetric_axis_frame(Z, atoms, tol, quantize_group_key=True)

    getcontext().prec -= 5

    return rot_mat


def orient_spherical(atoms, group, axes):
    """
    If all moments are equal, the system is of high symmetry. The principal axes of rotation are used.
    For tetrahedral symmetry, the 3 C2 axes are used as the Cartesian axes. For octahedral symmetry, the 3 C4
    axes are used. For icosahedral symmetry, an arbitrary C5 axis is chosen as the Z axis and the system is treated
    as a symmetric top.
    """
    rot_mat = SquareMatrix(3)
    tol = prec_tol()
    getcontext().prec += 2
    try:
        e_0 = Vector(3)
        e_1 = Vector(3)
        e_2 = Vector(3)
        e_0.assign(0, 1)
        e_1.assign(1, 1)
        e_2.assign(2, 1)

        #Td/Oh point groups: the three C2 (Td) or C4 (Oh) axes are used as the Cartesian axes
        if group in ('Td', 'Oh'):
            axis_z = max(axes, key = lambda a: abs(a.dot(e_2)))
            remaining = [a for a in axes if a is not axis_z]
            axis_y = max(remaining, key = lambda a: abs(a.dot(e_1)))
            axis_x = [a for a in remaining if a is not axis_y][0]

            rot_mat = rot_mat_from_axes([axis_x, axis_y, axis_z])

        #Ih point group: One C5 axis is chosen and used as Z. The system
        #is treated as a symmetric top.
        if group == 'Ih':
            Z = axes[0]
            rot_mat = symmetric_axis_frame(Z, atoms, tol, quantize_group_key=False)
    finally:
        getcontext().prec -= 2

    return rot_mat


def orient_asymm(eigvecs):
    """
    If all moments are unequal, the principal moments are used as the Cartesian axes.
    """
    return rot_mat_from_axes(eigvecs)


def standardize_axes(moments, eigvecs, atoms):
    """
    This function uses the principal moments of nuclear inertia to classify the type
    of top of the system: asymmetric, symmetric, spherical, linear, or a single
    atom. The system is then rotated to align the Cartesian axes with the principal
    axes.
    """

    tol = prec_tol()

    moment_a = moments[0]
    moment_b = moments[1]
    moment_c = moments[2]

    #Single atom
    if moment_a == moment_b == moment_c == 0:
        rot_mat = orient_atom()
    
    #Linear
    elif (moment_a == 0) and moment_b == moment_c:
        rot_mat = orient_linear(atoms)
    
    #Symmetric top
    elif (moment_a == moment_b != moment_c) or (moment_a != moment_b == moment_c):
        rot_mat = orient_symm(moment_a, moment_b, eigvecs, atoms)
    
    #Spherical top
    elif moment_a == moment_b == moment_c:
        group, axes = cn_axes_finder(atoms)
        rot_mat = orient_spherical(atoms, group, axes)
    
    #Asymmetric top
    elif moment_a != moment_b != moment_c:
        rot_mat = orient_asymm(eigvecs)

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

    if moment_a != moment_b != moment_c:
        standardized_atoms = fix_molecule_sign(standardized_atoms)
    getcontext().prec -= 10
    
    return standardized_atoms


def find_cn_axes(candidates, atoms, order, tol, extra_round=False):
    """
    Which of candidates are true order-fold rotation axes of atoms, tested
    via Rodrigues' rotation formula. extra_round matches cn_axes_finder's
    original C2 (mislabeled c3) test, which re-rounds v_rot one digit
    tighter before comparing; the C5/C4 tests never did this.
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

            if extra_round:
                rounded = Vector(3)
                getcontext().prec -= 1
                for i in range(3):
                    rounded.elements[i] = +v_rot.elements[i]
                getcontext().prec += 1
                v_rot = rounded

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


def cn_axes_finder(atoms):
    """
    This function is designed for the spherical top case. Finds principal rotation axes of high symmetry point groups.
    Identifies the system as tetrahedral, octahedral, or icosahedral and returns the principal axes.
    """

    getcontext().prec += 2
    try:
        tol = prec_tol(2)

        #Principal axes can only be either going through atoms or bisecting two atoms.
        #Find axes going through atoms
        candidates = []

        for atom in atoms:
            v = pos_vector(atom)
            norm = v.norm()
            if norm > tol:
                candidates.append(v.normalized())

        #Uniqueness check for axes through atoms
        unique = []
        for v in candidates:
            if not any(abs(abs(v.dot(u)) - 1) < tol for u in unique):
                unique.append(v)
        candidates = unique

        #Finds bisecting axes
        bisect_candidates = list(candidates)
        for v_0, v_1 in combinations(bisect_candidates, r = 2):
            w = v_0.add(v_1)
            candidates.append(w.normalized())

        #Another uniqueness check
        unique = []
        for v in candidates:
            if not any(abs(abs(v.dot(u)) - 1) < tol for u in unique):
                unique.append(v)
        candidates = unique

        #Rotation around the candidate principal axes using Rodrigues' rotation formula.
        #Icosahedral symmetry utilizes the C5 axes, octahedral C4, and tetrahedral C2 for their rotations.
        #Test for icosahedral:
        c5 = find_cn_axes(candidates, atoms, 5, tol)
        if len(c5) > 0:
            return 'Ih', c5

        #Test for octahedral:
        c4 = find_cn_axes(candidates, atoms, 4, tol)
        if len(c4) > 0:
            return 'Oh', c4

        #Test for tetrahedral:
        c2 = find_cn_axes(candidates, atoms, 2, tol, extra_round=True)
        if len(c2) > 0:
            return 'Td', c2
    finally:
        getcontext().prec -= 2


def atom_sort_key(atom):
    return (atom.charge, atom.x, atom.y, atom.z)


def fix_molecule_sign(atoms):

    choices = [
        (1, 1, 1),
        (1, -1, -1),
        (-1, 1, -1),
        (-1, -1, 1)
    ]

    best = None
    best_key = None

    for sx, sy, sz in choices:

        trial = []

        for atom in atoms:
            trial.append(
                Atom(
                    atom.element,
                    sx * atom.x,
                    sy * atom.y,
                    sz * atom.z,
                    atom.charge
                )
            )
        ordered = sorted(trial, key=atom_sort_key)

        key = tuple(atom_sort_key(atom) for atom in ordered)

        if best_key is None or key < best_key:
            best = ordered
            best_key = key

    return best