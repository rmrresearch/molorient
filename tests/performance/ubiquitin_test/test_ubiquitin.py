from molorient.cli import parse_xyz
from molorient.orient_system import orient_system


def test_ubiquitin(benchmark):
    atoms, folder, base, ext = parse_xyz("tests/performance/ubiquitin_test/1UBQ.xyz")
    oriented_atoms = orient_system(atoms)
    assert len(oriented_atoms) == len(atoms)
    benchmark(orient_system, atoms)
