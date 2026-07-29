"""NIPAM/DABCO copolymer — the two-arbitrary-residues generalization."""

from pathlib import Path

import pytest

from polybuilder import (
    PolymerSpec,
    build_polymer,
    load_user_library,
)
from polybuilder.exceptions import InvalidPolymerSpecError
from polybuilder.library import _reset_registry_for_tests

EXAMPLE = Path(__file__).resolve().parent.parent / "examples" / "cpp_monomers.py"


@pytest.fixture(autouse=True)
def _clean():
    _reset_registry_for_tests()
    load_user_library(str(EXAMPLE))
    yield
    _reset_registry_for_tests()


def test_nipam_first_dabco_comonomer_builds(tmp_path):
    """NIPAM as the majority residue, DABCO as the minority comonomer."""
    spec = PolymerSpec(
        first_residue="NIPAM",
        n_first=3,
        comonomer="DABCO",
        n_comonomer=1,
        repeats=2,
    )
    pdb = tmp_path / "p.pdb"
    rtp = tmp_path / "p.rtp"
    result = build_polymer(spec, str(pdb), str(rtp))

    # Sequence should be [NIPAM, NIPAM, NIPAM, DABCO] × 2
    assert result.sequence == ["NIPAM"] * 3 + ["DABCO"] + ["NIPAM"] * 3 + ["DABCO"]

    rtp_text = rtp.read_text()
    # NIPAM shows up as NIF/NIR/NIL, DABCO as DAR (or DAF/DAL depending on position)
    assert "[ NIF ]" in rtp_text
    assert "[ NIR ]" in rtp_text
    assert "[ DAR ]" in rtp_text or "[ DAL ]" in rtp_text
    # No MMA anywhere.
    assert "[ MCF ]" not in rtp_text
    assert "[ MCR ]" not in rtp_text


def test_nipam_dabco_with_kps_initiator(tmp_path):
    """The whole thermoresponsive+cationic recipe with sulfate end groups."""
    spec = PolymerSpec(
        first_residue="NIPAM",
        n_first=2,
        comonomer="DABCO",
        n_comonomer=1,
        repeats=2,
        cap=True,
        initiator="KPS",
        ends="both",
    )
    pdb = tmp_path / "p.pdb"
    rtp = tmp_path / "p.rtp"
    build_polymer(spec, str(pdb), str(rtp))
    rtp_text = rtp.read_text()
    assert "[ KPH ]" in rtp_text
    assert "[ KPT ]" in rtp_text
    assert "[ NIF ]" in rtp_text


def test_cap_uses_first_residue_not_mma(tmp_path):
    """--cap should append one first_residue unit, not hard-code MMA."""
    spec = PolymerSpec(
        first_residue="NIPAM",
        n_first=2,
        comonomer="DABCO",
        n_comonomer=1,
        repeats=1,
        cap=True,
    )
    pdb = tmp_path / "p.pdb"
    rtp = tmp_path / "p.rtp"
    result = build_polymer(spec, str(pdb), str(rtp))
    assert result.sequence[-1] == "NIPAM"
    assert "[ NIL ]" in rtp.read_text()   # last residue is NIPAM


def test_mma_defaults_still_work(tmp_path):
    """Not specifying first_residue keeps the legacy MMA-based behavior."""
    spec = PolymerSpec(
        comonomer="DMAPS", n_first=3, n_comonomer=1, repeats=2, extra_bridge=3,
    )
    assert spec.first_residue == "MMA"
    pdb = tmp_path / "p.pdb"
    rtp = tmp_path / "p.rtp"
    build_polymer(spec, str(pdb), str(rtp))
    text = rtp.read_text()
    assert "[ MCF ]" in text
    assert "[ DMR ]" in text


def test_legacy_n_mma_alias_still_works():
    """The old n_mma kwarg still populates n_first when first_residue is MMA."""
    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        spec = PolymerSpec(comonomer="DMAPS", n_mma=3, n_comonomer=1, repeats=1)
    assert spec.n_first == 3


def test_legacy_n_mma_rejected_with_non_mma_first_residue():
    """n_mma alias only applies when first_residue is MMA."""
    with pytest.raises(InvalidPolymerSpecError, match="n_mma"):
        PolymerSpec(
            first_residue="NIPAM",
            comonomer="DABCO",
            n_mma=3, n_comonomer=1, repeats=1,
        )


def test_bridge_rejected_when_no_sulfobetaine_present(tmp_path):
    """A pure NIPAM+DABCO polymer can't use extra_bridge > 0."""
    spec = PolymerSpec(
        first_residue="NIPAM", n_first=1,
        comonomer="DABCO", n_comonomer=1,
        repeats=1, extra_bridge=2,
    )
    with pytest.raises(InvalidPolymerSpecError, match="sulfobetaine"):
        build_polymer(spec, str(tmp_path / "p.pdb"), str(tmp_path / "p.rtp"))
