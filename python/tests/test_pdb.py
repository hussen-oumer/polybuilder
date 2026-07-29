from pathlib import Path

from polybuilder.pdb import interleave_pdb


def _write(path: Path, lines):
    path.write_text("\n".join(lines) + "\n")


def test_interleave_places_hydrogens_after_parent(tmp_path):
    src = tmp_path / "in.pdb"
    dst = tmp_path / "out.pdb"

    # Two heavy atoms in residue 1, each with one hydrogen written *after* the
    # heavies (not interleaved). interleave_pdb should regroup them.
    _write(src, [
        "ATOM      1 C1   MMC A   1       0.000   0.000   0.000  1.00  0.00           C",
        "ATOM      2 C2   MMC A   1       1.500   0.000   0.000  1.00  0.00           C",
        "ATOM      3 H1   MMC A   1       0.000   1.080   0.000  1.00  0.00           H",
        "ATOM      4 H2   MMC A   1       1.500   1.080   0.000  1.00  0.00           H",
    ])

    interleave_pdb(str(src), str(dst))

    atom_names = [
        line[12:16].strip()
        for line in dst.read_text().splitlines()
        if line.startswith("ATOM")
    ]
    assert atom_names == ["C1", "H1", "C2", "H2"]


def test_interleave_ignores_out_of_residue_hydrogens(tmp_path):
    src = tmp_path / "in.pdb"
    dst = tmp_path / "out.pdb"

    _write(src, [
        "ATOM      1 C1   MMC A   1       0.000   0.000   0.000  1.00  0.00           C",
        "ATOM      2 H1   MMC A   2       0.000   1.080   0.000  1.00  0.00           H",
    ])
    interleave_pdb(str(src), str(dst))
    body = dst.read_text().splitlines()
    atom_lines = [ln for ln in body if ln.startswith("ATOM")]
    # Only C1 is emitted; H1 is in a different residue and finds no parent.
    assert len(atom_lines) == 1
    assert atom_lines[0][12:16].strip() == "C1"


def test_interleave_treats_ch_named_carbons_as_heavy(tmp_path):
    """Regression: atom names like 'CH1' (aromatic C) must NOT be classified
    as hydrogens just because their second character is 'H'.  The element
    column (77-78) is the source of truth."""
    src = tmp_path / "in.pdb"
    dst = tmp_path / "out.pdb"

    _write(src, [
        "ATOM      1 CH1  RES A   1       0.000   0.000   0.000  1.00  0.00           C",
        "ATOM      2 CH2  RES A   1       1.400   0.000   0.000  1.00  0.00           C",
        "ATOM      3 CBZ  RES A   1       2.800   0.000   0.000  1.00  0.00           C",
        "ATOM      4 H1   RES A   1       0.000   1.080   0.000  1.00  0.00           H",
    ])
    interleave_pdb(str(src), str(dst))
    names = [
        ln[12:16].strip()
        for ln in dst.read_text().splitlines()
        if ln.startswith("ATOM")
    ]
    # All three C atoms must survive as heavy atoms.
    for c in ("CH1", "CH2", "CBZ"):
        assert c in names, f"heavy atom {c!r} was dropped by interleaver; got {names}"
