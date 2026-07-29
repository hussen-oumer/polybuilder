"""VBD / BVD residue tests (topology + CONECT records).

Counter-ions are not written by the builder — add them later during
solvation (Packmol / gmx genion).  These tests verify only the covalent
polymer topology.
"""

from pathlib import Path

import pytest

from polybuilder import (
    PolymerSpec,
    ResidueFragment,
    build_polymer,
    get_residue,
    load_user_library,
)
from polybuilder.config import AtomSpec
from polybuilder.library import _reset_registry_for_tests

EXAMPLE = Path(__file__).resolve().parent.parent / "examples" / "cpp_monomers.py"


@pytest.fixture(autouse=True)
def _clean():
    _reset_registry_for_tests()
    load_user_library(str(EXAMPLE))
    yield
    _reset_registry_for_tests()


def test_vbd_registered_as_residue_fragment():
    vbd = get_residue("VBD")
    assert isinstance(vbd, ResidueFragment)


def test_dabco_alias_points_to_same_residue():
    """The legacy DABCO name still resolves to VBD (backwards compatibility)."""
    assert get_residue("DABCO") is get_residue("VBD")


def test_vbd_uses_new_prefix_in_rtp(tmp_path):
    """A NIPAM + VBD copolymer emits VBx residue blocks, not DAx."""
    spec = PolymerSpec(
        first_residue="NIPAM", n_first=2,
        comonomer="VBD",      n_comonomer=1,
        repeats=2, cap=True,
    )
    pdb = tmp_path / "p.pdb"
    rtp = tmp_path / "p.rtp"
    build_polymer(spec, str(pdb), str(rtp))
    rtp_text = rtp.read_text()
    assert "[ VBR ]" in rtp_text or "[ VBF ]" in rtp_text or "[ VBL ]" in rtp_text
    assert "[ DAR ]" not in rtp_text


def test_dabco_legacy_name_uses_da_prefix(tmp_path):
    """Requesting comonomer='DABCO' keeps the historical DAx residue codes."""
    spec = PolymerSpec(
        first_residue="NIPAM", n_first=2,
        comonomer="DABCO",    n_comonomer=1,
        repeats=2, cap=True,
    )
    pdb = tmp_path / "p.pdb"
    rtp = tmp_path / "p.rtp"
    build_polymer(spec, str(pdb), str(rtp))
    rtp_text = rtp.read_text()
    assert "[ DAR ]" in rtp_text or "[ DAF ]" in rtp_text or "[ DAL ]" in rtp_text


def test_no_hetatm_ion_records_written(tmp_path):
    """No HETATM records for counter-ions should appear anywhere."""
    spec = PolymerSpec(
        first_residue="NIPAM", n_first=1,
        comonomer="VBD",      n_comonomer=1,
        repeats=3,
    )
    pdb = tmp_path / "p.pdb"
    rtp = tmp_path / "p.rtp"
    build_polymer(spec, str(pdb), str(rtp))
    assert "HETATM" not in pdb.read_text(), (
        "counter-ion HETATM records should NOT be written; add ions later "
        "with Packmol or gmx genion"
    )


def test_residuefragment_has_no_counter_ions_attribute():
    """After the removal, ResidueFragment should not expose counter_ions."""
    frag = ResidueFragment(
        smiles="C(C)",
        atoms=(AtomSpec("CC", "opls_140", 0.0, 1), AtomSpec("C1", "opls_135", 0.0, 1)),
    )
    assert not hasattr(frag, "counter_ions")


def test_pdb_contains_conect_records_for_cbz_n1_bond(tmp_path):
    """The critical vinyl-benzyl → DABCO bond (CBZ-N1) must appear as a
    CONECT record so viewers render it even if the C-N distance is a bit
    stretched by UFF."""
    spec = PolymerSpec(
        first_residue="MMA", n_first=1,
        comonomer="VBD",    n_comonomer=1,
        repeats=1,
    )
    pdb = tmp_path / "p.pdb"
    rtp = tmp_path / "p.rtp"
    build_polymer(spec, str(pdb), str(rtp))

    lines = pdb.read_text().splitlines()
    conect_lines = [ln for ln in lines if ln.startswith("CONECT")]
    assert conect_lines, "no CONECT records written"

    def serial_of(atom_name: str) -> int:
        for ln in lines:
            if ln.startswith("ATOM") and ln[12:16].strip() == atom_name:
                return int(ln[6:11])
        raise AssertionError(f"{atom_name} not found in PDB")

    cbz = serial_of("CBZ")
    n1 = serial_of("N1")

    def partners(central: int) -> set[int]:
        for ln in conect_lines:
            if int(ln[6:11]) == central:
                out = set()
                for start in (11, 16, 21, 26):
                    field = ln[start:start + 5]
                    if field.strip():
                        out.add(int(field))
                return out
        return set()

    assert n1 in partners(cbz), "CBZ CONECT missing N1"
    assert cbz in partners(n1), "N1 CONECT missing CBZ"


def test_bvd_has_two_benzenes_and_pendant_vinyl(tmp_path):
    """BVD (bis-benzyl-DABCO) has a second benzene ring on the pendant N⁺,
    with an unreacted vinyl at its para position."""
    spec = PolymerSpec(
        first_residue="MMA", n_first=1,
        comonomer="BVD",    n_comonomer=1,
        repeats=1,
    )
    pdb = tmp_path / "p.pdb"
    rtp = tmp_path / "p.rtp"
    build_polymer(spec, str(pdb), str(rtp))

    rtp_text = rtp.read_text()
    for atom in ("CG", "CD1", "CD2", "CE1", "CE2", "CZ"):
        assert f"  {atom}" in rtp_text, f"ring 1 atom {atom} missing"
    for atom in ("CG2", "CD3", "CD4", "CE3", "CE4", "CZ2"):
        assert f"  {atom}" in rtp_text, f"ring 2 atom {atom} missing"
    assert "  CBZ2   N2" in rtp_text or "  N2     CBZ2" in rtp_text
    assert "CV1" in rtp_text and "CV2" in rtp_text

    pdb_lines = pdb.read_text().splitlines()

    def serial_of(name: str) -> int:
        for ln in pdb_lines:
            if ln.startswith("ATOM") and ln[12:16].strip() == name:
                return int(ln[6:11])
        raise AssertionError(f"{name} not in PDB")

    def partners(central: int) -> set[int]:
        for ln in pdb_lines:
            if ln.startswith("CONECT") and int(ln[6:11]) == central:
                out: set[int] = set()
                for start in (11, 16, 21, 26):
                    field = ln[start:start + 5]
                    if field.strip():
                        out.add(int(field))
                return out
        return set()

    assert serial_of("CBZ2") in partners(serial_of("N2")), "N2-CBZ2 bond missing in CONECT"
    assert serial_of("CV2") in partners(serial_of("CV1")), "CV1-CV2 (pendant vinyl) missing"
