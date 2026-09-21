"""
Decimal-only eigensolvers for symmetric 3x3 matrices (Viete's trigonometric
method for eigenvalues, cross products for eigenvectors), plus the geometric
primitives shared by both orientation pipelines below.
"""

from decimal import Decimal, getcontext, ROUND_HALF_UP
from molorient.utils.trig_helpers import arccos_series, cos_series, pi_as_decimal
from molorient.utils.precision import prec_tol
from molorient.classes.vector import Vector
from molorient.classes.square_matrix import SquareMatrix


# Geometric primitives shared by both orientation pipelines (axis_standardization.py
# and nwchem_orientation.py both already import from this module, so hosting these
# here needs no new import cycle).
def pos_vector(atom):
    """An atom's (x, y, z) as a Vector(3)."""
    return Vector.from_atom(atom)


def build_inertia_tensor(atoms, weights):
    """
    Symmetric inertia tensor: sum_i weights[i] * (position outer-product
    formula). Shared by axis_standardization.inertia_tensor (charge-weighted)
    and nwchem_orientation.build_tensor (mass-weighted); weights supplies the
    per-atom scalar each one uses.
    """
    I_xx = I_yy = I_zz = I_xy = I_xz = I_yz = 0
    for w, a in zip(weights, atoms):
        I_xx += w * (a.y**2 + a.z**2)
        I_yy += w * (a.x**2 + a.z**2)
        I_zz += w * (a.x**2 + a.y**2)
        I_xy -= w * a.x * a.y
        I_xz -= w * a.x * a.z
        I_yz -= w * a.y * a.z

    tensor = SquareMatrix(3)
    tensor.elements = [[I_xx, I_xy, I_xz], [I_xy, I_yy, I_yz], [I_xz, I_yz, I_zz]]
    return tensor


def rot_mat_from_axes(axes):
    """axes[0], axes[1], axes[2] become the X, Y, Z columns of a SquareMatrix(3)."""
    return SquareMatrix.from_columns(axes)


def eigval_solver(squarematrix):
    """
    Finds eigenvalues of a 3x3 Hermitian Matrix using Viète's Trigonometric
    Method for the characteristic polynomial.
    """

    getcontext().prec += 5
    try:
        xx, xy, xz = squarematrix.elements[0]
        yy, yz = squarematrix.elements[1][1:]
        zz = squarematrix.elements[2][2]

        a = -1
        b = xx + yy + zz
        c = xy**2 + xz**2 + yz**2 - xx*yy - xx*zz - yy*zz
        d = xx*yy*zz + 2*xy*xz*yz - xx*yz**2 - yy*xz**2 - zz*xy**2

        # Converting to a depressed cubic form: t^3 + pt + q = 0
        # of which x = t - b/3a
        p = (3*a*c - b**2) / (3*a**2)
        q = (2*b**3 - 9*a*b*c + 27*a**2*d) / (27*a**3)

        tol = prec_tol(9)
        if abs(p) < tol:
            x_0 = x_1 = x_2 = b / 3

        else:
            z = (3 * q) / (2 * p) * (-3 / p).sqrt()
            z_rd = z.quantize(
                prec_tol(3), rounding = ROUND_HALF_UP
                )
            y = arccos_series(z_rd)

            t_0 = cos_series(y / 3)
            t_1 = cos_series((y / 3) - (2 * pi_as_decimal()) / 3)
            t_2 = cos_series((y / 3) - (4 * pi_as_decimal()) / 3)

            sqrt_term = (-p / 3).sqrt()
            coeff_term = -b / 3

            x_0 = (2 * sqrt_term * t_0) - coeff_term
            x_1 = (2 * sqrt_term * t_1) - coeff_term
            x_2 = (2 * sqrt_term * t_2) - coeff_term
    finally:
        getcontext().prec -= 5

    return +x_0, +x_1, +x_2


def best_candidate(candidates, eps):
    """
    The first of candidates whose squared length is within eps of the
    longest one, or None if the longest is zero (all candidates degenerate).
    """
    best_len = max(c.dot(c) for c in candidates)
    if best_len > 0:
        for c in candidates:
            if c.dot(c) >= eps * best_len:
                return c
    return None


def eigvec_solver(eig_0, eig_1, eig_2, squarematrix):
    """
    Solves for eigenvectors of 3x3 Hermitian matrix using the cross-product method to
    cut down on arithmetic operations.
    """

    getcontext().prec += 5
    try:
        e_0 = Vector(3)
        e_1 = Vector(3)
        e_2 = Vector(3)
        e_0.assign(0, Decimal('1.0'))
        e_1.assign(1, Decimal('1.0'))
        e_2.assign(2, 1)

        if eig_0 == eig_1 == eig_2 == Decimal(0):
            norm_vecs = [v.normalized() for v in [e_0, e_1, e_2]]

        else:
            a_0 = Vector(3)
            a_1 = Vector(3)
            a_2 = Vector(3)
            id_mat = SquareMatrix(3)
            cross_term_0 = Vector(3)
            cross_term_1 = Vector(3)
            cross_term_2 = Vector(3)

            for i in range(3):
                id_mat.elements[i][i] = Decimal('1.0')

            scaled_0 = e_0.scale(-eig_0)
            scaled_1 = e_1.scale(-eig_0)
            scaled_2 = e_2.scale(-eig_0)
            for i in range(3):
                a_0.elements[i] = squarematrix.elements[i][0]
                a_1.elements[i] = squarematrix.elements[i][1]
                a_2.elements[i] = squarematrix.elements[i][2]

                cross_term_0.elements[i] = a_0.elements[i] + scaled_0.elements[i]
                cross_term_1.elements[i] = a_1.elements[i] + scaled_1.elements[i]
                cross_term_2.elements[i] = a_2.elements[i] + scaled_2.elements[i]

            tol = prec_tol(2)
            eps = Decimal('1e-6')
            mus = []
            lin_ind = False

            for i in range(3):
                if cross_term_1.elements[i] != Decimal('0'):
                    mu = cross_term_0.elements[i] / cross_term_1.elements[i]
                    mus.append(mu)

                elif cross_term_0.elements[i] != Decimal('0'):
                    lin_ind = True
                    break

            if not lin_ind and (len(mus) != 0 and all(abs(m - mus[0]) < tol for m in mus)):
                scale_term = Decimal(1) / (Decimal(1) + mus[0]**2).sqrt()
                w = Vector(3)
                w.assign(0, 1)
                w.assign(1, -mus[0])
                v_0 = w.scale(scale_term)

            else:
                v0_candidates = [
                    cross_term_0.cross(cross_term_1),
                    cross_term_0.cross(cross_term_2),
                    cross_term_1.cross(cross_term_2),
                ]
                v_0 = best_candidate(v0_candidates, eps)
                if v_0 is None:
                    v_0 = Vector(3)
                    v_0.assign(0, 1)

            char_mat_0 = Vector(3)
            char_mat_1 = Vector(3)
            char_mat_2 = Vector(3)

            scaled_id = id_mat.scale(-eig_1)
            for i in range(3):
                char_mat_0.elements[i] = squarematrix.elements[i][0] + scaled_id.elements[i][0]
                char_mat_1.elements[i] = squarematrix.elements[i][1] + scaled_id.elements[i][1]
                char_mat_2.elements[i] = squarematrix.elements[i][2] + scaled_id.elements[i][2]

            v1_candidates = [
                v_0.cross(char_mat_0),
                v_0.cross(char_mat_1),
                v_0.cross(char_mat_2),
            ]
            v_1 = best_candidate(v1_candidates, eps)

            if v_1 is None:
                # v_0 != 0 here (the all-zero-eigenvalue case already
                # returned above), so at least one of v_0 x e_0/e_1/e_2 is
                # nonzero; take the largest instead of the first past tol,
                # which could otherwise settle on a near-zero vector.
                v_1 = max((v_0.cross(e) for e in [e_0, e_1, e_2]), key=lambda c: c.dot(c))

            v_2 = v_0.cross(v_1)

            norm_vecs = [v.normalized() for v in [v_0, v_1, v_2]]
    finally:
        getcontext().prec -= 5

    return norm_vecs