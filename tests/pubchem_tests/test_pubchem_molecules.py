from pubchem_functions.generate_molecules import search_pubchem, randomize_orientation
from molorient.utils.cli import parse_xyz, set_precision
from molorient.utils.orient_system import orient_system
from unittest.mock import patch


def test_pubchem_molecules_6_figs(tmp_path):
    count = 0
    with patch('builtins.input', return_value = '6'):
        result = set_precision()
        assert result == 6

    while count < 1000:
        atoms_output, cid = search_pubchem()
        count += 1
        filepath = tmp_path / f"{count}_pubchem.xyz"
        with open(filepath, 'w') as f:
            f.write(f"{len(atoms_output['atoms'])}\n")
            f.write(f"CID {cid} from PubChem\n")
            for atom in atoms_output['atoms']:
                f.write(f"{atom['element']} {atom['x']} {atom['y']} {atom['z']}\n")
        
        atoms, folder, base, ext = parse_xyz(filepath)
        std_atoms = orient_system(atoms)
        randomized_atoms = randomize_orientation(atoms)
        oriented_atoms = orient_system(randomized_atoms)

        for std_atom, oriented_atom in zip(std_atoms, oriented_atoms):
            assert std_atom.element == oriented_atom.element
            assert std_atom.x == oriented_atom.x
            assert std_atom.y == oriented_atom.y
            assert std_atom.z == oriented_atom.z

        count += 1


def test_pubchem_molecules_8_figs(tmp_path):
    count = 0
    with patch('builtins.input', return_value = '8'):
        result = set_precision()
        assert result == 8

    while count < 1000:
        atoms_output, cid = search_pubchem()
        count += 1
        filepath = tmp_path / f"{count}_pubchem.xyz"
        with open(filepath, 'w') as f:
            f.write(f"{len(atoms_output['atoms'])}\n")
            f.write(f"CID {cid} from PubChem\n")
            for atom in atoms_output['atoms']:
                f.write(f"{atom['element']} {atom['x']} {atom['y']} {atom['z']}\n")
        
        atoms, folder, base, ext = parse_xyz(filepath)
        std_atoms = orient_system(atoms)
        randomized_atoms = randomize_orientation(atoms)
        oriented_atoms = orient_system(randomized_atoms)

        for std_atom, oriented_atom in zip(std_atoms, oriented_atoms):
            assert std_atom.element == oriented_atom.element
            assert std_atom.x == oriented_atom.x
            assert std_atom.y == oriented_atom.y
            assert std_atom.z == oriented_atom.z

        count += 1


def test_pubchem_molecules_4_figs(tmp_path):
    count = 0
    with patch('builtins.input', return_value = '4'):
        result = set_precision()
        assert result == 4

    while count < 1000:
        atoms_output, cid = search_pubchem()
        count += 1
        filepath = tmp_path / f"{count}_pubchem.xyz"
        with open(filepath, 'w') as f:
            f.write(f"{len(atoms_output['atoms'])}\n")
            f.write(f"CID {cid} from PubChem\n")
            for atom in atoms_output['atoms']:
                f.write(f"{atom['element']} {atom['x']} {atom['y']} {atom['z']}\n")
        
        atoms, folder, base, ext = parse_xyz(filepath)
        for atom in atoms:
            atom.x = +atom.x
            atom.y = +atom.y
            atom.z = +atom.z

        std_atoms = orient_system(atoms)
        randomized_atoms = randomize_orientation(atoms)
        oriented_atoms = orient_system(randomized_atoms)

        for std_atom, oriented_atom in zip(std_atoms, oriented_atoms):
            assert std_atom.element == oriented_atom.element
            assert std_atom.x == oriented_atom.x
            assert std_atom.y == oriented_atom.y
            assert std_atom.z == oriented_atom.z

        count += 1