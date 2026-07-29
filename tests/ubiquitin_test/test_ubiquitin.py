from molorient.utils.cli import parse_xyz, set_precision
from molorient.utils.orient_system import orient_system


def test_ubiquitin(benchmark):
    atoms, folder, base, ext = parse_xyz("tests/ubiquitin_test/1UBQ.xyz")
    oriented_atoms = orient_system(atoms)
    benchmark(orient_system, atoms)
