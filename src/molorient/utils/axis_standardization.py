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

    #Calculate the principal moments of the inertia tensor
    eigvals, eigvecs = np.linalg.eig(inertia_tensor)

    #Sort the eigenvalues/ eigenvectorsin increasing order
    idx = np.argsort(eigvals)
    eigvals = eigvals[idx]
    eigvecs = eigvecs[:, idx]   
    norm_eigvecs = eigvecs / np.linalg.norm(eigvecs, axis=0)
    moment_a, moment_b, moment_c = eigvals

    return moment_a, moment_b, moment_c, norm_eigvecs


# def standardize_axes(atoms, moment_a, moment_b, moment_c, norm_eigvecs):
#     """
#     This function standardizes the axes of the system by rotating the coordinates of the atoms.
#     The axes are defined by whether the system is an asymmetric top, symmetric top, spherical top, linear or a 
#     single atom based on the principal moments. The function returns a new list of Atom objects with the standardized coordinates.
#     """

#     #Asymmetric top: all three moments are different
#     if moment_a != moment_b and moment_b != moment_c and moment_a != moment_c:
#         rotation_matrix = norm_eigvecs
    
#     #Single atom: all three moments are zero
#     elif moment_a == 0 and moment_b == 0 and moment_c == 0:
#         rotation_matrix = np.eye(3)
    
#     #Linear molecule: one moment is zero and the other two are equal
#     elif (moment_a == 0 and moment_b == moment_c) or (moment_b == 0 and moment_a == moment_c) or (moment_c == 0 and moment_a == moment_b):
#         vec = np.array([
#             atoms[-1].x - atoms[0].x,
#             atoms[-1].y - atoms[0].y,
#             atoms[-1].z - atoms[0].z
#         ])

#         vec /= np.linalg.norm(vec)
#         theta = np.arccos(vec[2])
#         phi = np.arctan2(vec[1], vec[0])

#         rot_matrix = np.array([
#             [np.cos(theta), -np.sin(phi), 0],
#             [np.cos(theta) * np.sin(phi), np.cos(theta) * np.cos(phi), -np.sin(phi)],
#             [np.sin(theta) * np.sin(phi), np.sin(theta) * np.cos(phi), np.cos(theta)]
#         ])
#         rotation_matrix = rot_matrix.T
    
#     #Spherical top: all three moments are equal

#     #Symmetric top: two moments are equal and the third is different
#     #This case evaluates every atom based on distance to XY plane, positive Z projection,
#     #distance to Z axis, and the lowest atomic number.
#     else:
#         if moment_a != moment_b:
#             Z = norm_eigvecs[:, 2]
#         elif moment_a != moment_c:
#             Z = norm_eigvecs[:, 1]
#         else:
#             Z = norm_eigvecs[:, 0]
    
#         Z = Z / np.linalg.norm(Z)

#         r = np.array([
#             [atom.x for atom in atoms],
#             [atom.y for atom in atoms],
#             [atom.z for atom in atoms]
#         ])

#         # Build circular set groups
#         groups = {}
#         for atom in atoms:
#             pos = [atom.x, atom.y, atom.z]
#             z_proj = np.dot(Z, pos)
#             dist_to_z = np.linalg.norm(pos - z_proj * Z)
#             key = round(float(z_proj), 8)
#             if key not in groups:
#                 groups[key] = {
#                     'atoms': [],
#                     'z_proj': z_proj,
#                     'dist_to_xy': abs(z_proj),
#                     'dist_to_z': dist_to_z,
#                 }
#             groups[key]['atoms'].append(atom)

#             candidates = list(groups.values())
            
#             #1. Nearest to xy plane
#             min_dist_to_xy = min(g['dist_to_xy'] for g in candidates)
#             candidates = [g for g in candidates if np.isclose(g['dist_to_xy'], min_dist_to_xy, atol=1e-8)]

#             #2. Positive Z proj
#             pos_z = [g for g in candidates if g['z_proj'] > 0]
#             if pos_z:
#                 candidates = pos_z
            
#             #3. Nearest to Z axis
#             min_dist_to_z = min(g['dist_to_z'] for g in candidates)
#             candidates = [g for g in candidates if np.isclose(g['dist_to_z'], min_dist_to_z, atol=1e-8)]

#             #4. Lowest atomic number
#             min_charge = min(min(atom.charge for atom in g['atoms']) for g in candidates)
#             candidates = [g for g in candidates if any(np.isclose(atom.charge, min_charge) for atom in g['atoms'])]

#             key_set = candidates[0]

#             key_atom = key_set['atoms'][0]
#             key_pos = np.array([key_atom.x, key_atom.y, key_atom.z])
#             key_proj = np.dot(Z, key_pos)
#             Y = key_pos - key_proj * Z
#             Y = Y / np.linalg.norm(Y)
#             X = np.cross(Y, Z)
#             rotation_matrix = np.column_stack((X, Y, Z))

#     standardized_atoms = []

#     for atom in atoms:
#         pos = np.array([atom.x, atom.y, atom.z])
#         new_pos = rotation_matrix.T @ pos
#         standardized_atoms.append(Atom(atom.element, new_pos[0], new_pos[1], new_pos[2], atom.mass, atom.charge))

#     return standardized_atoms
