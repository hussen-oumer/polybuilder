import pytest

from polybuilder.exceptions import InvalidSmilesError
from polybuilder.pubchem import generate_config_entry


def test_methacrylate_gets_alpha_methyl_type():
    smiles, atoms = generate_config_entry("MMA", "C=C(C)C(=O)OC")
    assert smiles is not None
    assert atoms is not None
    assert atoms[0].name == "CC"
    assert atoms[0].ff_type == "opls_139"  # α-methyl vinyl


def test_acrylate_gets_non_methyl_type():
    smiles, atoms = generate_config_entry("acryl", "C=CC(=O)OC")
    assert smiles is not None
    assert atoms[0].ff_type == "opls_140"


def test_non_polymerisable_returns_none():
    smiles, atoms = generate_config_entry("EtOH", "CCO")
    assert smiles is None
    assert atoms is None


def test_invalid_smiles_raises():
    with pytest.raises(InvalidSmilesError):
        generate_config_entry("bad", "this is not smiles !!!")

    with pytest.raises(InvalidSmilesError):
        generate_config_entry("empty", "")
