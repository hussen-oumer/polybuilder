"""End-to-end builds of the single-fragment monomers from cpp_monomers.py."""

from pathlib import Path

import pytest

from polybuilder import PolymerSpec, build_polymer, load_user_library
from polybuilder.exceptions import InvalidPolymerSpecError
from polybuilder.library import _reset_registry_for_tests

EXAMPLE = Path(__file__).resolve().parent.parent / "examples" / "cpp_monomers.py"


@pytest.fixture(autouse=True)
def _reset_registry():
    _reset_registry_for_tests()
    load_user_library(str(EXAMPLE))
    yield
    _reset_registry_for_tests()


@pytest.mark.parametrize("comonomer", ["NIPAM", "IBOA", "DABCO"])
def test_builds_copolymer_with_single_fragment_comonomer(tmp_path, comonomer):
    """A small MMA/comonomer copolymer builds and emits both files."""
    spec = PolymerSpec(
        comonomer=comonomer,
        n_mma=1,
        n_comonomer=1,
        repeats=2,
        extra_bridge=0,
    )
    pdb = tmp_path / f"{comonomer}.pdb"
    rtp = tmp_path / f"{comonomer}.rtp"
    result = build_polymer(spec, str(pdb), str(rtp))

    assert pdb.stat().st_size > 0
    assert rtp.stat().st_size > 0
    assert result.sequence == ["MMA", comonomer, "MMA", comonomer]
    # Should contain a residue block for the new comonomer prefix.
    rtp_text = rtp.read_text()
    assert "[ MMC ]" in rtp_text
    assert f"[ {comonomer[:2]}F ]" in rtp_text or f"[ {comonomer[:2]}R ]" in rtp_text


def test_builds_nipam_homopolymer(tmp_path):
    """No MMA at all — pure NIPAM chain."""
    spec = PolymerSpec(
        comonomer="NIPAM", n_mma=0, n_comonomer=4, repeats=1,
    )
    pdb = tmp_path / "n.pdb"
    rtp = tmp_path / "n.rtp"
    build_polymer(spec, str(pdb), str(rtp))
    text = rtp.read_text()
    assert "[ NIF ]" in text
    assert "[ NIR ]" in text
    assert "[ NIL ]" in text


def test_bridge_rejected_for_single_fragment(tmp_path):
    """extra_bridge > 0 with a single-fragment comonomer is a spec error."""
    spec = PolymerSpec(
        comonomer="NIPAM", n_mma=1, n_comonomer=1, repeats=1, extra_bridge=2,
    )
    with pytest.raises(InvalidPolymerSpecError, match="sulfobetaine"):
        build_polymer(spec, str(tmp_path / "n.pdb"), str(tmp_path / "n.rtp"))
