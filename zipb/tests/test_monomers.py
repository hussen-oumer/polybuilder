import pytest

from polybuilder.support.exceptions import UnknownMonomerError
from polybuilder.library.monomers import (
    build_dynamic_a3316,
    build_dynamic_dmaps,
    build_dynamic_residue,
)


def test_dmaps_matches_generic_lookup():
    smiles_a, atoms_a = build_dynamic_dmaps(3)
    smiles_b, atoms_b = build_dynamic_residue("DMAPS", 3)
    assert smiles_a == smiles_b
    assert [a.as_tuple() for a in atoms_a] == [a.as_tuple() for a in atoms_b]


def test_a3316_bridge_count_matches_expected_atoms():
    _, atoms_zero = build_dynamic_a3316(0)
    _, atoms_five = build_dynamic_a3316(5)
    assert len(atoms_five) == len(atoms_zero) + 5


def test_unknown_residue_raises():
    with pytest.raises(UnknownMonomerError):
        build_dynamic_residue("NOPE", 3)


def test_single_fragment_residue_requires_zero_bridge():
    with pytest.raises(ValueError):
        build_dynamic_residue("MMA", 3)
