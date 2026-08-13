"""Tests for the simple homopolymer() / copolymer() API helpers."""

import pytest

from polybuilder import build_polymer, copolymer, homopolymer
from polybuilder.support.exceptions import InvalidPolymerSpecError
from polybuilder.library import _reset_registry_for_tests


@pytest.fixture(autouse=True)
def _clean():
    _reset_registry_for_tests()
    yield
    _reset_registry_for_tests()


def test_homopolymer_spec():
    spec = homopolymer("DMAPS", 5)
    assert spec.first_residue == "DMAPS"
    assert spec.n_first == 5
    assert spec.n_comonomer == 0
    assert spec.repeats == 1


def test_homopolymer_builds(tmp_path):
    spec = homopolymer("DMAPS", 4)
    result = build_polymer(
        spec, str(tmp_path / "p.pdb"), str(tmp_path / "p.rtp")
    )
    assert result.sequence == ["DMAPS"] * 4


def test_homopolymer_initiator_and_bridge():
    spec = homopolymer("DMAPS", 3, initiator="KPS", ends="head", bridge=2)
    assert spec.initiator == "KPS"
    assert spec.ends == "head"
    assert spec.extra_bridge == 2


def test_copolymer_alternating():
    spec = copolymer("DMAPS", "A3367", 3)
    assert (spec.n_first, spec.n_comonomer, spec.repeats) == (1, 1, 3)


def test_copolymer_block():
    spec = copolymer("DMAPS", "A3367", 3, block=True)
    assert (spec.n_first, spec.n_comonomer, spec.repeats) == (3, 3, 1)


def test_copolymer_explicit_overrides_n():
    spec = copolymer("DMAPS", "A3367", 99, na=3, nb=1, repeats=4)
    assert (spec.n_first, spec.n_comonomer, spec.repeats) == (3, 1, 4)


def test_copolymer_alternating_sequence(tmp_path):
    spec = copolymer("DMAPS", "A3367", 2)
    result = build_polymer(
        spec, str(tmp_path / "p.pdb"), str(tmp_path / "p.rtp")
    )
    assert result.sequence == ["DMAPS", "A3367", "DMAPS", "A3367"]


def test_copolymer_requires_counts():
    with pytest.raises(InvalidPolymerSpecError):
        copolymer("DMAPS", "A3367")
