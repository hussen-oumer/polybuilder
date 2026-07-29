"""Runtime registry of monomer residues.

The built-in residues in :mod:`polybuilder.config` (MMA, DMAPS, A3316) are
registered at import time.  Additional residues can be added at runtime
without touching the package:

* Programmatically::

      from polybuilder import register_residue, BetaineBase, ResidueFragment
      register_residue("MY_MONO", ResidueFragment(...))            # single fragment
      register_residue("MY_BET",  BetaineBase(head=..., tail=...)) # betaine

* From a user file (Python module)::

      polybuilder.load_user_library("~/monomers/mine.py")

  The file must expose either a ``MONOMERS`` dict (``name`` → residue) or
  call :func:`register_residue` at import time.  It is executed in a fresh
  module namespace so it never mutates the package on disk.

Legacy names ``register_betaine`` / ``get_betaine`` / ``available_betaines``
remain as aliases so pre-0.3 user files still load.
"""

from __future__ import annotations

import importlib.util
import logging
import os
import warnings
from collections.abc import Iterable
from pathlib import Path

from .config import (
    INITIATOR_LIBRARY,
    RESIDUE_LIBRARY,
    BetaineBase,
    Initiator,
    Residue,
    ResidueFragment,
)
from .exceptions import PolybuilderError, UnknownMonomerError

log = logging.getLogger(__name__)


# Runtime registries, seeded with the built-in libraries.  User plugins mutate
# copies of these dicts at load time.
_REGISTRY: dict[str, Residue] = dict(RESIDUE_LIBRARY)
_INITIATOR_REGISTRY: dict[str, Initiator] = dict(INITIATOR_LIBRARY)


def register_residue(name: str, residue: Residue, *, overwrite: bool = False) -> None:
    """Add or replace a residue in the runtime registry.

    ``residue`` may be a :class:`ResidueFragment` (single-fragment monomer
    such as NIPAM or DABCO) or a :class:`BetaineBase` (head + tunable CH2
    bridge + tail).

    Raises :class:`PolybuilderError` if ``name`` already exists and
    ``overwrite`` is ``False``.
    """
    if not isinstance(residue, (ResidueFragment, BetaineBase)):
        raise TypeError(
            "register_residue expected a ResidueFragment or BetaineBase, "
            f"got {type(residue).__name__}"
        )
    if name in _REGISTRY and not overwrite:
        raise PolybuilderError(
            f"residue '{name}' already registered; pass overwrite=True to replace."
        )
    _REGISTRY[name] = residue
    log.info("Registered residue '%s' (%d built-in, %d total)",
             name, len(RESIDUE_LIBRARY), len(_REGISTRY))


def get_residue(name: str) -> Residue:
    """Look up a residue by name; raise :class:`UnknownMonomerError` if absent."""
    try:
        return _REGISTRY[name]
    except KeyError as exc:
        known = ", ".join(sorted(_REGISTRY)) or "<none>"
        raise UnknownMonomerError(
            f"Unknown residue '{name}'. Known: {known}."
        ) from exc


def available_residues() -> Iterable[str]:
    """Return the currently registered residue names (sorted)."""
    return sorted(_REGISTRY)


# ---------------------------------------------------------------------------
# Initiators
# ---------------------------------------------------------------------------
def register_initiator(name: str, initiator: Initiator, *, overwrite: bool = False) -> None:
    """Add or replace an initiator in the runtime registry."""
    if not isinstance(initiator, Initiator):
        raise TypeError(
            f"register_initiator expected an Initiator, got {type(initiator).__name__}"
        )
    if name in _INITIATOR_REGISTRY and not overwrite:
        raise PolybuilderError(
            f"initiator '{name}' already registered; pass overwrite=True to replace."
        )
    _INITIATOR_REGISTRY[name] = initiator
    log.info(
        "Registered initiator '%s' (%d built-in, %d total)",
        name, len(INITIATOR_LIBRARY), len(_INITIATOR_REGISTRY),
    )


def get_initiator(name: str) -> Initiator:
    """Look up an initiator by name; raise :class:`UnknownMonomerError` if absent."""
    try:
        return _INITIATOR_REGISTRY[name]
    except KeyError as exc:
        known = ", ".join(sorted(_INITIATOR_REGISTRY)) or "<none>"
        raise UnknownMonomerError(
            f"Unknown initiator '{name}'. Known: {known}."
        ) from exc


def available_initiators() -> Iterable[str]:
    """Return the currently registered initiator names (sorted)."""
    return sorted(_INITIATOR_REGISTRY)


# ---------------------------------------------------------------------------
# Deprecated aliases kept for backward compatibility with 0.2 user files.
# ---------------------------------------------------------------------------
def register_betaine(name: str, base: Residue, *, overwrite: bool = False) -> None:
    warnings.warn(
        "register_betaine is deprecated; use register_residue.",
        DeprecationWarning,
        stacklevel=2,
    )
    register_residue(name, base, overwrite=overwrite)


def get_betaine(name: str) -> Residue:
    warnings.warn(
        "get_betaine is deprecated; use get_residue.",
        DeprecationWarning,
        stacklevel=2,
    )
    return get_residue(name)


def available_betaines() -> Iterable[str]:
    warnings.warn(
        "available_betaines is deprecated; use available_residues.",
        DeprecationWarning,
        stacklevel=2,
    )
    return available_residues()


def load_user_library(path: os.PathLike | str, *, overwrite: bool = False) -> Iterable[str]:
    """Load additional monomers and/or initiators from a user Python file.

    The file may:

    1. Expose ``MONOMERS = {"NAME": <ResidueFragment or BetaineBase>, ...}``,
    2. Expose ``INITIATORS = {"NAME": Initiator(...), ...}``, and/or
    3. Call :func:`register_residue` / :func:`register_initiator` at module
       top level (also the deprecated :func:`register_betaine`).

    All styles can coexist.  Returns the combined list of names newly
    registered by this load, prefixed with ``residue:`` or ``initiator:``
    to disambiguate.
    """
    file = Path(path).expanduser().resolve()
    if not file.is_file():
        raise FileNotFoundError(f"user library not found: {file}")

    spec = importlib.util.spec_from_file_location(f"polybuilder_user_{file.stem}", file)
    if spec is None or spec.loader is None:
        raise PolybuilderError(f"could not load user library: {file}")

    module = importlib.util.module_from_spec(spec)
    before_res = set(_REGISTRY)
    before_ini = set(_INITIATOR_REGISTRY)
    try:
        spec.loader.exec_module(module)
    except Exception as exc:
        raise PolybuilderError(
            f"error while executing user library {file}: {exc}"
        ) from exc

    monomers = getattr(module, "MONOMERS", None)
    if isinstance(monomers, dict):
        for name, residue in monomers.items():
            register_residue(name, residue, overwrite=overwrite)

    initiators = getattr(module, "INITIATORS", None)
    if isinstance(initiators, dict):
        for name, init in initiators.items():
            register_initiator(name, init, overwrite=overwrite)

    added_res = sorted(set(_REGISTRY) - before_res)
    added_ini = sorted(set(_INITIATOR_REGISTRY) - before_ini)
    added = [f"residue:{n}" for n in added_res] + [f"initiator:{n}" for n in added_ini]

    if not added:
        log.warning("user library %s did not register any new monomers or initiators", file)
    else:
        log.info("Loaded from %s: %s", file, ", ".join(added))
    return added


def _reset_registry_for_tests() -> None:
    """Restore the registries to the built-in libraries only. Intended for tests."""
    _REGISTRY.clear()
    _REGISTRY.update(RESIDUE_LIBRARY)
    _INITIATOR_REGISTRY.clear()
    _INITIATOR_REGISTRY.update(INITIATOR_LIBRARY)


__all__ = [
    "register_residue",
    "get_residue",
    "available_residues",
    "register_initiator",
    "get_initiator",
    "available_initiators",
    "load_user_library",
    # deprecated aliases
    "register_betaine",
    "get_betaine",
    "available_betaines",
]
