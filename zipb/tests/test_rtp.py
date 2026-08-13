
from polybuilder.write.rtp import reorder_rtp_by_cgnr


def test_reorder_rtp_by_cgnr_sorts_atoms(tmp_path):
    src = tmp_path / "in.rtp"
    dst = tmp_path / "out.rtp"

    src.write_text(
        "[ bondedtypes ]\n"
        "; junc  ang  dih  imp\n"
        " 1 1 1 1\n"
        "\n"
        "[ MMC ]\n"
        " [ atoms ]\n"
        "  A1   opls_1     0.0000   3\n"
        "  A2   opls_1     0.0000   1\n"
        "  A3   opls_1     0.0000   2\n"
        "\n [ bonds ]\n"
        "  A1     A2\n"
        "\n"
    )
    reorder_rtp_by_cgnr(str(src), str(dst))
    text = dst.read_text()

    # Atoms should appear sorted by their 4th column (cgnr).  Restrict to
    # lines with 4 fields to skip the bond line "A1     A2".
    atom_lines = [
        ln.strip() for ln in text.splitlines()
        if ln.strip().startswith(("A1", "A2", "A3")) and len(ln.split()) >= 4
    ]
    cgnrs = [int(ln.split()[3]) for ln in atom_lines]
    assert cgnrs == [1, 2, 3]
    # Bond section is preserved after atoms.
    assert "[ bonds ]" in text
    assert "A1     A2" in text
