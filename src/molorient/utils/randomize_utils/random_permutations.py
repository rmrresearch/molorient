from molorient.utils.orient_system import orient_system
import random


def random_permutation(atoms):
    shuffled_atoms = random.sample(atoms, len(atoms))

    return shuffled_atoms