from molorient.utils.cli import parse_xyz, set_precision
from unittest.mock import patch


def test_read_xyz(tmp_path):
    test_file = tmp_path / "test_xyz.xyz"

    test_file.write_text(
        "4\n"
        "Comment line\n"
        "H 0 0 0\n"
        "H 1 0 0\n"
        "H 0 0 1\n"
        "H 0 1 0\n"
    )

    atoms, folder, base, ext = parse_xyz(test_file)
    assert len(atoms) == 4
    assert atoms[0].element == 'H'
    assert atoms[0].x == 0
    assert folder == str(tmp_path)
    assert base == 'test_xyz'
    assert ext == '.xyz'


def test_get_precision():
    with patch('builtins.input', return_value = '10'):
        result = set_precision()
        assert result == 10