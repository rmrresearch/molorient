from decimal import Decimal, getcontext, ROUND_HALF_UP
from itertools import combinations
from molorient.classes.atom import Atom
from molorient.classes.square_matrix import SquareMatrix
from molorient.classes.vector import Vector
from molorient.utils.diagonalization import eigval_solver, eigvec_solver
from molorient.utils.trig_helpers import sin_series, cos_series, arccos_series, arctan2, pi_as_decimal


def inertia_tensor(atoms):
    """
    This function calculates the inertia tensor of the system using charge instead of mass. It returns
    the principal moments of inertia (eigenvalues) of the tensor in order of increasing magnitude.
    """

    tensor = SquareMatrix(3)
    I_xx = sum([atom.charge * (atom.y**2 + atom.z**2) for atom in atoms])
    I_yy = sum([atom.charge * (atom.x**2 + atom.z**2) for atom in atoms])
    I_zz = sum([atom.charge * (atom.x**2 + atom.y**2) for atom in atoms])
    I_xy = -sum([atom.charge * atom.x * atom.y for atom in atoms])
    I_xz = -sum([atom.charge * atom.x * atom.z for atom in atoms])
    I_yz = -sum([atom.charge * atom.y * atom.z for atom in atoms])

    tensor.assign(0, 0, I_xx)
    tensor.assign(0, 1, I_xy)
    tensor.assign(0, 2, I_xz)
    tensor.assign(1, 0, I_xy)
    tensor.assign(1, 1, I_yy)
    tensor.assign(1, 2, I_yz)
    tensor.assign(2, 0, I_xz)
    tensor.assign(2, 1, I_yz)
    tensor.assign(2, 2, I_zz)

    moment_a, moment_b, moment_c = sorted(eigval_solver(tensor))
    v_0, v_1, v_2 = eigvec_solver(moment_a, moment_b, moment_c, tensor)
    eigvals = [moment_a, moment_b, moment_c]
    eigvecs = [v_0, v_1, v_2]
    
    #Round values
    #(1) If number of decimal places > number of sig figs, round to 10**(-precision)
    #(2) If number of sig figs < number of decimal places, round to narest precision
    for i, e in enumerate(eigvals):
        t = e.as_tuple()
        sig_figs = len(t.digits)
        dec_places = max(0, -t.exponent)
        if sig_figs < dec_places:
            eigvals[i] = e.quantize(Decimal(10)**-(getcontext().prec))
        else:
            eigvals[i] = +e

    return eigvals, eigvecs


def standardize_axes(moments, eigvecs, atoms):
    """
    This function uses the principal moments of nuclear inertia to classify the type
    of top of the system: asymmetric, symmetric, spherical, linear, or a single
    atom. The system is then rotated to align the Cartesian axes with the principal
    axes.
    """

    tol = Decimal(10)**-(getcontext().prec)

    moment_a = moments[0]
    moment_b = moments[1]
    moment_c = moments[2]
    rot_mat = SquareMatrix(3)

    #Single atom: all three moments are zero
    if all(m == 0 for m in moments):
        for i in range(3):
            rot_mat.elements[i][i] = Decimal('1')

    #Linear: one moment is 0 and the others are equal
    elif (moment_a == 0) and moment_b == moment_c:
        getcontext().prec += 2
        
        vec = Vector(3)
        vec.elements[0] = (atoms[-1].x - atoms[0].x)
        vec.elements[1] = (atoms[-1].y - atoms[0].y)
        vec.elements[2] = (atoms[-1].z - atoms[0].z)

        norm = Decimal('1') / (vec.elements[0]**2 + vec.elements[1]**2 + vec.elements[2]**2).sqrt()
        norm_vec = vec.scale(norm)

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
    
    #Symmetric top: two moments are equal and the third is different.
    #This case evaluates every group of atoms based on its distance to the XY plane,
    #projection and distance to the new z axis, and the lowest atomic number.
    elif (moment_a == moment_b != moment_c) or (moment_a != moment_b == moment_c):
        getcontext().prec += 2
        tol = Decimal(10)**-(getcontext().prec - 2)

        if moment_a == moment_b:
            z_col = 2
        else:
            z_col = 0
        Z = Vector(3)
        for i in range(3):
            Z.elements[i] = eigvecs[z_col].elements[i]

        groups = {}
        for atom in atoms:
            pos = Vector(3)
            pos.elements[0] = atom.x
            pos.elements[1] = atom.y
            pos.elements[2] = atom.z
            Z_proj = Z.dot(pos)
            along = Z.scale(Z_proj)
            perp = Vector(3)
            for i in range(3):
                perp.elements[i] = pos.elements[i] - along.elements[i]
            dist_to_z = (perp.dot(perp)).sqrt()
            if dist_to_z == 0:
                continue

            key = (atom.element, Z_proj)
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
        key_pos = Vector(3)
        key_pos.elements[0] = key_atom.x
        key_pos.elements[1] = key_atom.y
        key_pos.elements[2] = key_atom.z

        key_Z_proj = Z.dot(key_pos)
        Y = key_pos.add(Z.scale(-key_Z_proj))
        norm = Decimal('1') / (Y.elements[0]**2 + Y.elements[1]**2 + Y.elements[2]**2).sqrt()
        Y_norm = Y.scale(norm)
        X = Y_norm.cross(Z)
        X_norm = Decimal('1') / (X.elements[0]**2 + X.elements[1]**2 + X.elements[2]**2).sqrt()
        X_norm = X.scale(X_norm)

        for i in range(3):
            rot_mat.elements[i][0] = X_norm.elements[i]
            rot_mat.elements[i][1] = Y_norm.elements[i]
            rot_mat.elements[i][2] = Z.elements[i]
            rot_mat.transpose()
        
        getcontext().prec -= 2
    
    #Spherical top: all three moments are the same.
    elif moment_a == moment_b == moment_c:
        group, axes = cn_axes_finder(atoms)
        getcontext().prec += 2
        
        e_0 = Vector(3)
        e_1 = Vector(3)
        e_2 = Vector(3)
        e_0.assign(0, 1)
        e_1.assign(1, 1)
        e_2.assign(2, 1)

        #Td point group: the three C2 axes are used as the Cartesian axes
        if group == 'Td':
            c2_z = max(axes, key = lambda a: abs(a.dot(e_2)))
            remaining = [a for a in axes if a is not c2_z]
            c2_y = max(remaining, key = lambda a: abs(a.dot(e_1)))
            c2_x = [a for a in remaining if a is not c2_y][0]

            for i in range(3):
                rot_mat.elements[i][0] = c2_x.elements[i]
                rot_mat.elements[i][1] = c2_y.elements[i]
                rot_mat.elements[i][2] = c2_z.elements[i]
            
            getcontext().prec -= 2
        
        #Oh point group: the three C4 axes are used as the Cartesian axes
        if group == 'Oh':
            c4_z = max(axes, key = lambda a: abs(a.dot(e_2)))
            remaining = [a for a in axes if a is not c4_z]
            c4_y = max(remaining, key = lambda a: abs(a.dot(e_1)))
            c4_x = [a for a in remaining if a is not c4_y][0]

            for i in range(3):
                rot_mat.elements[i][0] = c4_x.elements[i]
                rot_mat.elements[i][1] = c4_y.elements[i]
                rot_mat.elements[i][2] = c4_z.elements[i]

            getcontext().prec -= 2
            
        #Ih point group: One C5 axis is chosen and used as Z. The system
        #is treated as a symmetric top.
        if group == 'Ih':
            axes[0] = Z
            groups = {}
            for atom in atoms:
                pos = Vector(3)
                pos.elements[0] = atom.x
                pos.elements[1] = atom.y
                pos.elements[2] = atom.z
                Z_proj = Z.dot(pos)
                along = Z.scale(Z_proj)
                perp = Vector(3)
                for i in range(3):
                    perp.elements[i] = pos.elements[i] - along.elements[i]
                dist_to_z = (perp.dot(perp)).sqrt()
                if dist_to_z == 0:
                    continue

                key = (atom.element, Z_proj)
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
            key_pos = Vector(3)
            key_pos.elements[0] = key_atom.x
            key_pos.elements[1] = key_atom.y
            key_pos.elements[2] = key_atom.z

            key_Z_proj = Z.dot(key_pos)
            Y = key_pos.add(Z.scale(-key_Z_proj))
            norm = Decimal('1') / (Y.elements[0]**2 + Y.elements[1]**2 + Y.elements[2]**2).sqrt()
            Y_norm = Y.scale(norm)
            X = Y_norm.cross(Z)
            X_norm = Decimal('1') / (X.elements[0]**2 + X.elements[1]**2 + X.elements[2]**2).sqrt()
            X_norm = X.scale(X_norm)

            for i in range(3):
                rot_mat.elements[i][0] = X_norm.elements[i]
                rot_mat.elements[i][1] = Y_norm.elements[i]
                rot_mat.elements[i][2] = Z.elements[i]
                rot_mat.transpose()
                
    #Asymmetric top: all three moments are different
    elif moment_a != moment_b != moment_c:
        for i in range(3):
                rot_mat.elements[i][0] = eigvecs[0].elements[i]
                rot_mat.elements[i][1] = eigvecs[1].elements[i]
                rot_mat.elements[i][2] = eigvecs[2].elements[i]
    
    #Rotation
    standardized_atoms = []
    getcontext().prec += 1
    tol = Decimal(1).scaleb(-(getcontext().prec - 1))

    for atom in atoms:
        pos_vec = Vector(3)
        pos_vec.assign(0, atom.x)
        pos_vec.assign(1, atom.y) 
        pos_vec.assign(2, atom.z)

        new_pos = (rot_mat.transpose()).multiply(pos_vec)

        standardized_atoms.append(Atom(atom.element, 
                                        new_pos.elements[0].quantize(tol, rounding = ROUND_HALF_UP),
                                        new_pos.elements[1].quantize(tol, rounding = ROUND_HALF_UP),
                                        new_pos.elements[2].quantize(tol, rounding = ROUND_HALF_UP),
                                        atom.mass, 
                                        atom.charge))

    getcontext().prec -= 1

    return standardized_atoms


def cn_axes_finder(atoms):
    """
    This function is designed for the spherical top case. Finds principal rotation axes of high symmetry point groups.
    Identifies the system as tetrahedral, octahedral, or icosahedral and returns the principal axes.
    """

    getcontext().prec += 2
    tol = Decimal(10)**-(getcontext().prec - 2)

    #Principal axes can only be either going through atoms or bisecting two atoms.
    #Find axes going through atoms
    candidates = []

    for atom in atoms:
        v = Vector(3)
        v.elements[0] = atom.x
        v.elements[1] = atom.y
        v.elements[2] = atom.z
        norm = (v.dot(v)).sqrt()
        if norm > tol:
            candidates.append(v.scale(1 / norm))

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
        norm = (w.dot(w)).sqrt()
        candidates.append(w.scale(1 / norm))

    #Another uniqueness check
    unique = []
    for v in candidates:
        if not any(abs(abs(v.dot(u)) - 1) < tol for u in unique):
            unique.append(v)
    candidates = unique

    #Rotation around the candidate principal axes using Rodrigues' rotation formula.
    #Icosahedral symmetry utilizes the C5 axes, octahedral C4, and tetrahedral C2 for their rotations.
    #Test for icosahedral:
    theta = 2 * pi_as_decimal() / 5
    cos_theta = cos_series(theta)
    sin_theta = sin_series(theta)
    c5 = []
    for k in candidates:
        is_c5 = True
        for atom in atoms:
            v = Vector(3)
            v.elements[0] = atom.x
            v.elements[1] = atom.y
            v.elements[2] = atom.z
            v_rot = v.scale(cos_theta).add((k.cross(v)).scale(sin_theta).add(k.scale((k.dot(v)) * (1 - cos_theta))))
            
            if not any(
                (v_rot.elements[0] - atom2.x)**2 +
                (v_rot.elements[1] - atom2.y)**2 +
                (v_rot.elements[2] - atom2.z)**2 < tol
                for atom2 in atoms
            ):
                is_c5 = False
                break

        if is_c5:
            getcontext().prec -= 2
            c5.append(k)
    
    if len(c5) > 0:
        getcontext().prec -= 2
        return 'Ih', c5
    
    #Test for octahedral:
    theta = 2 * pi_as_decimal() / 4
    cos_theta = cos_series(theta)
    sin_theta = sin_series(theta)
    c4 = []
    for k in candidates:
        is_c4 = True
        for atom in atoms:
            v = Vector(3)
            v.elements[0] = atom.x
            v.elements[1] = atom.y
            v.elements[2] = atom.z
            v_rot = v.scale(cos_theta).add((k.cross(v)).scale(sin_theta).add(k.scale((k.dot(v)) * (1 - cos_theta))))
            
            if not any(
                (v_rot.elements[0] - atom2.x)**2 +
                (v_rot.elements[1] - atom2.y)**2 +
                (v_rot.elements[2] - atom2.z)**2 < tol
                for atom2 in atoms
            ):
                is_c4 = False
                break

        if is_c4:
            c4.append(k)

    if len(c4) > 0:
        getcontext().prec -= 2
        return 'Oh', c4

    #Test for tetrahedral:
    theta = 2 * pi_as_decimal() / 2
    cos_theta = cos_series(theta)
    sin_theta = sin_series(theta)
    c3 = []
    for k in candidates:
        is_c3 = True
        for atom in atoms:
            v = Vector(3)
            v.elements[0] = atom.x
            v.elements[1] = atom.y
            v.elements[2] = atom.z
            v_rot = v.scale(cos_theta).add((k.cross(v)).scale(sin_theta).add(k.scale((k.dot(v)) * (1 - cos_theta))))
            v_rot_rd = Vector(3)
            getcontext().prec -= 1
            for i in range(3):
                v_rot_rd.elements[i] = +v_rot.elements[i]
            getcontext().prec += 1

            if not any(
                (v_rot_rd.elements[0] - atom2.x)**2 +
                (v_rot_rd.elements[1] - atom2.y)**2 +
                (v_rot_rd.elements[2] - atom2.z)**2 < tol
                for atom2 in atoms
            ):
                is_c3 = False
                break

        if is_c3:
            c3.append(k)
    
    if len(c3) > 0:
        for v in c3:
            getcontext().prec -= 2
            return 'Td', c3