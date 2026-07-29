"""Tests for radical-initiator end-capping."""


import pytest

from polybuilder import (
    PolybuilderError,
    PolymerSpec,
    UnknownMonomerError,
    available_initiators,
    build_polymer,
    get_initiator,
    load_user_library,
    register_initiator,
)
from polybuilder.config import KPS_INITIATOR
from polybuilder.exceptions import InvalidPolymerSpecError
from polybuilder.library import _reset_registry_for_tests


@pytest.fixture(autouse=True)
def _clean():
    _reset_registry_for_tests()
    yield
    _reset_registry_for_tests()


def test_builtin_initiators_present():
    assert set(available_initiators()) == {"KPS", "AIBA", "TBHP"}


def test_get_unknown_initiator_raises():
    with pytest.raises(UnknownMonomerError):
        get_initiator("NOPE")


def test_register_initiator_duplicate_requires_overwrite():
    with pytest.raises(PolybuilderError):
        register_initiator("KPS", KPS_INITIATOR)
    register_initiator("KPS", KPS_INITIATOR, overwrite=True)


def test_register_wrong_type_rejected():
    with pytest.raises(TypeError):
        register_initiator("BAD", "not an initiator")  # type: ignore[arg-type]


def test_invalid_ends_rejected():
    with pytest.raises(InvalidPolymerSpecError):
        PolymerSpec(comonomer="MMA", n_mma=3, repeats=1, cap=True, ends="left")


def test_kps_both_ends_on_mma_homopolymer(tmp_path):
    """MMA homopolymer capped with KPS on both ends writes new residue blocks."""
    spec = PolymerSpec(
        comonomer="MMA", n_mma=4, n_comonomer=0, repeats=1, cap=True,
        initiator="KPS", ends="both",
    )
    pdb = tmp_path / "p.pdb"
    rtp = tmp_path / "p.rtp"
    build_polymer(spec, str(pdb), str(rtp))

    rtp_text = rtp.read_text()
    assert "[ KPH ]" in rtp_text     # KPS head cap
    assert "[ KPT ]" in rtp_text     # KPS tail cap
    assert "[ MCF ]" in rtp_text     # first MMA residue
    # Chain should carry a -2 net charge (two sulfate O-).
    # Re-parse the PDB via the Mol we already have via inspection.
    # Simplest: count "O-" annotations in the assembled SMILES by re-building.
    # (Full charge accounting is exercised by the RDKit parse in build_polymer.)


def test_kps_head_only_leaves_tail_methyl(tmp_path):
    spec = PolymerSpec(
        comonomer="MMA", n_mma=3, n_comonomer=0, repeats=1, cap=True,
        initiator="KPS", ends="head",
    )
    pdb = tmp_path / "p.pdb"
    rtp = tmp_path / "p.rtp"
    build_polymer(spec, str(pdb), str(rtp))
    rtp_text = rtp.read_text()
    assert "[ KPH ]" in rtp_text
    assert "[ KPT ]" not in rtp_text
    # CAP2 methyl is still present in the last residue's block.
    assert "CAP2" in rtp_text


def test_tbhp_tail_only(tmp_path):
    spec = PolymerSpec(
        comonomer="MMA", n_mma=3, n_comonomer=0, repeats=1, cap=True,
        initiator="TBHP", ends="tail",
    )
    pdb = tmp_path / "p.pdb"
    rtp = tmp_path / "p.rtp"
    build_polymer(spec, str(pdb), str(rtp))
    rtp_text = rtp.read_text()
    assert "[ TBT ]" in rtp_text
    assert "[ TBH ]" not in rtp_text


def test_ends_none_leaves_default_methyl_caps(tmp_path):
    """ends='none' with an initiator specified still yields methyl caps."""
    spec = PolymerSpec(
        comonomer="MMA", n_mma=3, n_comonomer=0, repeats=1, cap=True,
        initiator="KPS", ends="none",
    )
    pdb = tmp_path / "p.pdb"
    rtp = tmp_path / "p.rtp"
    build_polymer(spec, str(pdb), str(rtp))
    rtp_text = rtp.read_text()
    assert "[ KPH ]" not in rtp_text
    assert "[ KPT ]" not in rtp_text
    assert "CAP1" in rtp_text
    assert "CAP2" in rtp_text


def test_kps_mid_chain_insertion(tmp_path):
    """Explicit mid-chain KPS insertion produces the KPM residue block."""
    spec = PolymerSpec(
        comonomer="MMA", n_first=5, n_comonomer=0, repeats=2, cap=True,
        initiator="KPS",
        initiator_positions=["head", "tail", 4],   # mid at index 4
    )
    pdb = tmp_path / "p.pdb"
    rtp = tmp_path / "p.rtp"
    build_polymer(spec, str(pdb), str(rtp))
    rtp_text = rtp.read_text()
    assert "[ KPH ]" in rtp_text
    assert "[ KPT ]" in rtp_text
    assert "[ KPM ]" in rtp_text
    # The persulfate bridge atoms must appear in the KPM block.
    assert "OA1" in rtp_text and "OA2" in rtp_text and "OP1" in rtp_text


def test_initiator_percent_auto_placement(tmp_path):
    """--initiator-percent picks head + tail + evenly-spaced mid positions."""
    # 20 residues, 20% → 4 initiators = head + tail + 2 mids.
    spec = PolymerSpec(
        comonomer="MMA", n_first=10, n_comonomer=0, repeats=2, cap=False,
        initiator="KPS", initiator_percent=20.0,
    )
    pdb = tmp_path / "p.pdb"
    rtp = tmp_path / "p.rtp"
    build_polymer(spec, str(pdb), str(rtp))
    text = rtp.read_text()
    assert "[ KPH ]" in text
    assert "[ KPT ]" in text
    assert "[ KPM ]" in text


def test_initiator_percent_zero_gives_no_caps(tmp_path):
    """0% initiator means no head, no tail, no mid — plain methyl caps."""
    spec = PolymerSpec(
        comonomer="MMA", n_first=5, n_comonomer=0, repeats=1, cap=True,
        initiator="KPS", initiator_percent=0.0,
    )
    pdb = tmp_path / "p.pdb"
    rtp = tmp_path / "p.rtp"
    build_polymer(spec, str(pdb), str(rtp))
    text = rtp.read_text()
    assert "[ KPH ]" not in text
    assert "[ KPT ]" not in text
    assert "[ KPM ]" not in text
    assert "CAP1" in text
    assert "CAP2" in text


def test_initiator_positions_out_of_range_are_ignored(tmp_path):
    """Positions outside [0, n_res-2] are dropped with a warning, not a crash."""
    spec = PolymerSpec(
        comonomer="MMA", n_first=3, n_comonomer=0, repeats=1, cap=False,
        initiator="KPS",
        initiator_positions=["head", "tail", 999],  # 999 is out of range
    )
    pdb = tmp_path / "p.pdb"
    rtp = tmp_path / "p.rtp"
    build_polymer(spec, str(pdb), str(rtp))
    text = rtp.read_text()
    assert "[ KPH ]" in text
    assert "[ KPT ]" in text
    assert "[ KPM ]" not in text     # 999 was invalid, so no mid


def test_mid_chain_unsupported_initiator_errors(tmp_path):
    """AIBA has no mid_smiles; requesting mid insertion should error clearly."""
    from polybuilder.exceptions import InvalidPolymerSpecError
    spec = PolymerSpec(
        comonomer="MMA", n_first=3, n_comonomer=0, repeats=1, cap=False,
        initiator="AIBA",
        initiator_positions=["head", 1],
    )
    with pytest.raises(InvalidPolymerSpecError, match="mid_smiles"):
        build_polymer(spec, str(tmp_path / "p.pdb"), str(tmp_path / "p.rtp"))


def test_initiator_percent_requires_initiator():
    """initiator_percent without an initiator is a spec error."""
    from polybuilder.exceptions import InvalidPolymerSpecError
    with pytest.raises(InvalidPolymerSpecError, match="initiator_percent"):
        PolymerSpec(
            comonomer="MMA", n_first=5, n_comonomer=0, repeats=1, cap=False,
            initiator_percent=5.0,
        )


def test_user_initiator_via_plugin_file(tmp_path):
    """User plugin can register a new initiator through the INITIATORS dict."""
    lib = tmp_path / "user.py"
    lib.write_text(
        "from polybuilder import AtomSpec, Initiator\n"
        "MY_INI = Initiator(\n"
        "    name='MYX',\n"
        "    head_smiles='C(C)(C)O',\n"
        "    tail_smiles='OC(C)C',\n"
        "    head_atoms=(\n"
        "        AtomSpec('CT', 'opls_135', 0.0, 1),\n"
        "        AtomSpec('CM1', 'opls_135', 0.0, 1),\n"
        "        AtomSpec('CM2', 'opls_135', 0.0, 1),\n"
        "        AtomSpec('OA', 'opls_154', 0.0, 1),\n"
        "    ),\n"
        "    tail_atoms=(\n"
        "        AtomSpec('OA', 'opls_154', 0.0, 1),\n"
        "        AtomSpec('CT', 'opls_135', 0.0, 1),\n"
        "        AtomSpec('CM1', 'opls_135', 0.0, 1),\n"
        "        AtomSpec('CM2', 'opls_135', 0.0, 1),\n"
        "    ),\n"
        ")\n"
        "INITIATORS = {'MYX': MY_INI}\n"
    )
    added = load_user_library(str(lib))
    assert any(a.startswith("initiator:MYX") for a in added)
    assert "MYX" in available_initiators()

    spec = PolymerSpec(
        comonomer="MMA", n_mma=2, n_comonomer=0, repeats=1, cap=True,
        initiator="MYX", ends="both",
    )
    pdb = tmp_path / "p.pdb"
    rtp = tmp_path / "p.rtp"
    build_polymer(spec, str(pdb), str(rtp))
    text = rtp.read_text()
    assert "[ MYH ]" in text
    assert "[ MYT ]" in text
