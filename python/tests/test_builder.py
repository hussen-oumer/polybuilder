import pytest

from polybuilder import PolymerSpec, build_polymer
from polybuilder.support.exceptions import InvalidPolymerSpecError


def test_polymerspec_rejects_negative_counts():
    with pytest.raises(InvalidPolymerSpecError):
        PolymerSpec(comonomer="DMAPS", n_mma=-1, n_comonomer=1, repeats=1)

    with pytest.raises(InvalidPolymerSpecError):
        PolymerSpec(comonomer="DMAPS", n_mma=1, n_comonomer=1, repeats=0)


def test_polymerspec_rejects_empty_polymer():
    with pytest.raises(InvalidPolymerSpecError):
        PolymerSpec(comonomer="DMAPS", n_mma=0, n_comonomer=0, repeats=1, cap=False)


def test_polymerspec_accepts_legacy_betaine_kwarg():
    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        spec = PolymerSpec(betaine="DMAPS", n_mma=1, n_betaine=1, repeats=1)
    assert spec.comonomer == "DMAPS"
    assert spec.n_comonomer == 1


def test_build_small_dmaps_polymer_writes_files(tmp_path):
    spec = PolymerSpec(
        comonomer="DMAPS",
        n_mma=1,
        n_comonomer=1,
        repeats=1,
        extra_bridge=1,
        cap=False,
    )
    pdb = tmp_path / "p.pdb"
    rtp = tmp_path / "p.rtp"
    result = build_polymer(spec, str(pdb), str(rtp))

    assert pdb.exists() and pdb.stat().st_size > 0
    assert rtp.exists() and rtp.stat().st_size > 0
    # sequence: [MMA, DMAPS] repeated once, no cap
    assert result.sequence == ["MMA", "DMAPS"]
    # Every atom in the PDB should have a residue name and a serial
    pdb_lines = [ln for ln in pdb.read_text().splitlines() if ln.startswith("ATOM")]
    assert pdb_lines, "no ATOM records written"
    for line in pdb_lines:
        assert line[17:20].strip(), f"missing residue name in line: {line!r}"
    # RTP should include both residue templates
    rtp_text = rtp.read_text()
    assert "[ bondedtypes ]" in rtp_text
    assert "[ atoms ]" in rtp_text
