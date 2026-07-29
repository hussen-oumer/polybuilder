from pathlib import Path

import pytest

from polybuilder import (
    PolybuilderError,
    UnknownMonomerError,
    available_residues,
    build_polymer,
    get_residue,
    load_user_library,
    register_residue,
)
from polybuilder.config import DMAPS_BASE
from polybuilder.library import _reset_registry_for_tests


@pytest.fixture(autouse=True)
def _clean_registry():
    _reset_registry_for_tests()
    yield
    _reset_registry_for_tests()


def test_builtins_present():
    # MMA now lives in the same registry as the sulfobetaines so the builder
    # can look it up polymorphically.
    assert set(available_residues()) == {"MMA", "DMAPS", "A3316"}


def test_register_and_get_new_betaine():
    register_residue("MYMONO", DMAPS_BASE)
    assert "MYMONO" in available_residues()
    assert get_residue("MYMONO") is DMAPS_BASE


def test_register_wrong_type_rejected():
    with pytest.raises(TypeError):
        register_residue("BAD", object())  # type: ignore[arg-type]


def test_double_register_requires_overwrite():
    register_residue("MYMONO", DMAPS_BASE)
    with pytest.raises(PolybuilderError):
        register_residue("MYMONO", DMAPS_BASE)
    # explicit overwrite is allowed
    register_residue("MYMONO", DMAPS_BASE, overwrite=True)


def test_unknown_lookup_raises():
    with pytest.raises(UnknownMonomerError):
        get_residue("NOPE")


def test_load_user_library_dict_style(tmp_path):
    lib = tmp_path / "userlib.py"
    lib.write_text(
        "from polybuilder.config import DMAPS_BASE\n"
        "MONOMERS = {'MYCOPY': DMAPS_BASE}\n"
    )
    added = load_user_library(str(lib))
    assert list(added) == ["residue:MYCOPY"]
    assert get_residue("MYCOPY") is DMAPS_BASE


def test_load_user_library_register_style(tmp_path):
    lib = tmp_path / "userlib.py"
    lib.write_text(
        "from polybuilder import register_residue\n"
        "from polybuilder.config import DMAPS_BASE\n"
        "register_residue('MYCOPY', DMAPS_BASE)\n"
    )
    added = load_user_library(str(lib))
    assert "residue:MYCOPY" in added


def test_load_missing_file():
    with pytest.raises(FileNotFoundError):
        load_user_library("/does/not/exist/anywhere.py")


def test_end_to_end_with_plugin_monomer(tmp_path):
    """A build using a plugin-registered monomer succeeds and writes files."""
    repo_examples = Path(__file__).resolve().parent.parent / "examples" / "cpp_monomers.py"
    assert repo_examples.is_file(), "example plugin file missing"
    load_user_library(str(repo_examples))

    from polybuilder import PolymerSpec

    pdb = tmp_path / "p.pdb"
    rtp = tmp_path / "p.rtp"
    build_polymer(
        PolymerSpec(
            comonomer="A3367", n_mma=1, n_comonomer=1, repeats=1, extra_bridge=1,
        ),
        str(pdb), str(rtp),
    )
    assert pdb.stat().st_size > 0
    assert rtp.stat().st_size > 0
