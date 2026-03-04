from geom_standardizer.classes.atom import Atom
from geom_standardizer.utils.translation import translation_vector, translate_to_origin

def test_translation():
    atoms = [
        Atom("H", 0.0, 0.0, 0.0, 1.008, 1.0),
        Atom("H", 1.0, 0.0, 0.0, 1.008, 1.0),
    ]

    expected_coords = [
        (-0.5, 0.0, 0.0),
        (0.5, 0.0, 0.0),
    ]
    
    trans_vector = translation_vector(atoms)
    translated_atoms = translate_to_origin(atoms, trans_vector)

    for atom, (x_exp, y_exp, z_exp) in zip(translated_atoms, expected_coords):
        assert atom.x == x_exp
        assert atom.y == y_exp
        assert atom.z == z_exp