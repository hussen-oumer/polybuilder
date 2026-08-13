"""Monomer & initiator registry (config, catalog, runtime registration).

Importing this package registers the built-in residues and initiators.
The public registry API is re-exported here so callers can use
``from polybuilder.library import get_residue`` without knowing that the
implementation lives in :mod:`polybuilder.library.registry`.
"""

from __future__ import annotations

from .registry import (
    _reset_registry_for_tests,
    available_betaines,
    available_initiators,
    available_residues,
    get_betaine,
    get_initiator,
    get_residue,
    load_user_library,
    register_betaine,
    register_initiator,
    register_residue,
)

__all__ = [
    "register_residue",
    "get_residue",
    "available_residues",
    "register_initiator",
    "get_initiator",
    "available_initiators",
    "load_user_library",
    "register_betaine",
    "get_betaine",
    "available_betaines",
]
