from polybuilder.config import (
    A3316_BASE,
    BETAINE_LIBRARY,
    DMAPS_BASE,
    MMA_RESIDUE,
    AtomSpec,
)


def test_library_contains_known_betaines():
    assert set(BETAINE_LIBRARY) == {"DMAPS", "A3316"}


def test_mma_atoms_have_expected_names():
    names = [a.name for a in MMA_RESIDUE.atoms]
    assert names == ["CC", "CB", "CD", "OD2", "OD1", "CE"]


def test_dmaps_bridge_extends_smiles_and_atoms():
    base = DMAPS_BASE
    smiles0, atoms0 = base.with_bridge(0)
    smiles2, atoms2 = base.with_bridge(2)

    assert smiles0 == base.head.smiles + base.tail.smiles
    assert smiles2 == base.head.smiles + "CC" + base.tail.smiles
    assert len(atoms2) == len(atoms0) + 2
    bridge_atoms = [a for a in atoms2 if a.name.startswith("CX")]
    assert [a.name for a in bridge_atoms] == ["CX1", "CX2"]
    for a in bridge_atoms:
        assert isinstance(a, AtomSpec)
        assert a.ff_type == "opls_136"


def test_a3316_bridge_cgnr_matches_base():
    _, atoms = A3316_BASE.with_bridge(1)
    bridge = next(a for a in atoms if a.name == "CX1")
    assert bridge.cgnr == 6


def test_negative_bridge_rejected():
    import pytest
    with pytest.raises(ValueError):
        DMAPS_BASE.with_bridge(-1)
