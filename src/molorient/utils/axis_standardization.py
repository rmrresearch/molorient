from decimal import Decimal
from molorient.classes.atom import Atom
import numpy as np


def inertia_tensor(atoms):
    """
    This function calculates the inertia tensor of the system using charge instead of mass. It returns
    the principal moments of inertia (eigenvalues) of the tensor in order of increasing magnitude.
    """

    #Inertia tensor
    I_xx = sum(float(atom.charge * (atom.y**2 + atom.z**2)) for atom in atoms)
    I_yy = sum(float(atom.charge * (atom.x**2 + atom.z**2)) for atom in atoms)
    I_zz = sum(float(atom.charge * (atom.x**2 + atom.y**2)) for atom in atoms)
    I_xy = -sum(float(atom.charge * atom.x * atom.y) for atom in atoms)
    I_xz = -sum(float(atom.charge * atom.x * atom.z) for atom in atoms)
    I_yz = -sum(float(atom.charge * atom.y * atom.z) for atom in atoms)

    inertia_tensor = np.array([[I_xx, I_xy, I_xz],
                               [I_xy, I_yy, I_yz],
                               [I_xz, I_yz, I_zz]])

    # Calculate the principal moments of the inertia tensor
    eigvals, eigvecs = np.linalg.eigh(inertia_tensor)

    # Sort the eigenvalues/eigenvectors in increasing order
    idx = np.argsort(eigvals)
    eigvals = eigvals[idx]
    eigvecs = eigvecs[:, idx]

    # Normalize eigenvectors
    norm_eigvecs = eigvecs / np.linalg.norm(eigvecs, axis=0)

    # Force largest component of each eigenvector to be positive
    for i in range(norm_eigvecs.shape[1]):
        idx_max = np.argmax(np.abs(norm_eigvecs[:, i]))
        if norm_eigvecs[idx_max, i] < 0:
            norm_eigvecs[:, i] *= -1

    moment_a, moment_b, moment_c = eigvals

    return moment_a, moment_b, moment_c, norm_eigvecs


def standardize_axes(atoms, moment_a, moment_b, moment_c, norm_eigvecs):
    """
    This function standardizes the axes of the system by rotating the coordinates of the atoms.
    The axes are defined by whether the system is an asymmetric top, symmetric top, spherical top, linear or a 
    single atom based on the principal moments. The function returns a new list of Atom objects with the standardized coordinates.
    """
        
    print("Principal moments of inertia:")
    print(f"  A: {moment_a}")
    print(f"  B: {moment_b}")
    print(f"  C: {moment_c}")

    print("Normalized eigenvectors (columns):")
    print(norm_eigvecs)

    #Single atom: all three moments are zero
    if moment_a == 0 and moment_b == 0 and moment_c == 0:
        rotation_matrix = np.eye(3)

    #Asymmetric top: all three moments are different
    elif np.isclose(moment_a, moment_b, atol=1e-6) == False and np.isclose(moment_b, moment_c, atol=1e-6) == False and np.isclose(moment_a, moment_c, atol=1e-6) == False:
        rotation_matrix = norm_eigvecs
    
    #Linear molecule: one moment is zero and the other two are equal
    elif (np.isclose(moment_a, 0) and np.isclose(moment_b, moment_c)):
        vec = np.array([
            float(atoms[-1].x - atoms[0].x),
            float(atoms[-1].y - atoms[0].y),
            float(atoms[-1].z - atoms[0].z)
        ])
        vec /= np.linalg.norm(vec)
        
        theta = np.arccos(vec[2])
        phi = np.arctan2(vec[1], vec[0])

        rotation_matrix = np.array([
            [0.5*np.cos(theta + phi) + 0.5*np.cos(theta - phi), np.sin(phi), 0.5*np.sin(theta + phi) + 0.5*np.sin(theta - phi)],
            [-0.5*np.sin(theta + phi) + 0.5*np.sin(theta - phi), np.cos(phi), 0.5*np.cos(theta + phi) - 0.5*np.cos(theta - phi)],
            [-np.sin(theta), 0, np.cos(theta)]
        ])

#     #Symmetric top: two moments are equal and the third is different
#     #This case evaluates every atom based on distance to XY plane, positive Z projection,
#     #distance to Z axis, and the low008est atomic number.
    elif (np.isclose(moment_a, moment_b, atol=1e-6) and not np.isclose(moment_b, moment_c, atol=1e-6)) or (np.isclose(moment_a, moment_c, atol=1e-6) and not np.isclose(moment_c, moment_b, atol=1e-6)) or (np.isclose(moment_b, moment_c, atol=1e-6) and not np.isclose(moment_c, moment_a, atol=1e-6)):
        try:
            if not np.isclose(moment_a, moment_b, atol=1e-6):
                Z = norm_eigvecs[:, 1]
            elif not np.isclose(moment_a, moment_c, atol=1e-6):
                Z = norm_eigvecs[:, 2]
            else:
                Z = norm_eigvecs[:, 0]

        except Exception as e:
            print("Error determining Z axis for symmetric top:", e)
        # Build circular set groups
        groups = {}
        for atom in atoms:
            pos = [float(atom.x), float(atom.y), float(atom.z)]
            z_proj = np.dot(Z, pos)
            dist_to_z = np.linalg.norm(pos - z_proj * Z)
            if np.isclose(dist_to_z, 0, atol=1e-8):
                continue
            key = (atom.element, round(float(z_proj), 8))
            if key not in groups:
                groups[key] = {
                    'atoms': [],
                    'z_proj': z_proj,
                    'dist_to_xy': abs(z_proj),
                    'dist_to_z': dist_to_z,
                }
            groups[key]['atoms'].append(atom)

        candidates = list(groups.values())
        
        #1. Nearest to xy plane
        min_dist_to_xy = min(g['dist_to_xy'] for g in candidates)
        candidates = [g for g in candidates if np.isclose(g['dist_to_xy'], min_dist_to_xy, atol=1e-8)]

        #2. Positive Z proj
        pos_z = [g for g in candidates if g['z_proj'] > 0]
        if pos_z:
            candidates = pos_z
        
        #3. Nearest to Z axis
        min_dist_to_z = min(g['dist_to_z'] for g in candidates)
        candidates = [g for g in candidates if np.isclose(g['dist_to_z'], min_dist_to_z, atol=1e-8)]

        #4. Lowest atomic number
        min_charge = min(min(float(atom.charge) for atom in g['atoms']) for g in candidates)
        candidates = [g for g in candidates if any(np.isclose(float(atom.charge), float(min_charge)) for atom in g['atoms'])]

        key_set = candidates[0]

        key_atom = key_set['atoms'][0]
        key_pos = np.array([float(key_atom.x), float(key_atom.y), float(key_atom.z)])
        key_proj = np.dot(Z, key_pos)
        Y = key_pos - key_proj * Z
        Y = Y / np.linalg.norm(Y)
        X = np.cross(Y, Z)
        rotation_matrix = np.column_stack((X, Y, Z))
        print("Key atom for symmetric top:", key_pos)
        print("Rotation matrix for symmetric top:", rotation_matrix)

#     #Spherical top: all three moments are equal

    standardized_atoms = []

    for atom in atoms:
        pos = np.array([float(atom.x), float(atom.y), float(atom.z)])
        new_pos = rotation_matrix.T @ pos
        decimal_pos = [Decimal(str(coord)) for coord in new_pos]
        standardized_atoms.append(Atom(atom.element, decimal_pos[0], decimal_pos[1], decimal_pos[2], atom.mass, atom.charge))
    
    print("Rotation matrix:")
    print(rotation_matrix)

    print("Standardized coordinates:")
    for atom in standardized_atoms:
        print(f"  {atom.element}: ({atom.x}, {atom.y}, {atom.z})")

    return standardized_atoms
