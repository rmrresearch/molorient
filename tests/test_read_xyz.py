from molorient.utils.read_xyz import read_xyz


def test_read_xyz(tmp_path):
    test_file = tmp_path / "test_xyz.xyz"

    test_file.write_text(
        "4"
        "\n Comment line"
        "\nH 0 0 0"
        "\nH 1 0 0"
        "\nH 0 0 1"
        "\nH 0 1 0"
    )

    atoms = read_xyz(test_file)
    assert len(atoms) == 4
    assert atoms[0].element == 'H'
    assert atoms[0].x == 0
    assert 2 == 1