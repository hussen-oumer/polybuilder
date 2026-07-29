"""Uniform residue-lookup entry point.

Whether a residue is a sulfobetaine (head + tunable CH2 bridge + tail) or a
single-fragment monomer (MMA, NIPAM, DABCO, ...), :func:`build_dynamic_residue`
returns a ``(smiles, atoms)`` pair suitable for the polymer builder.

Legacy helper names ``build_dynamic_dmaps``, ``build_dynamic_a3316`` and
``build_dynamic_betaine`` remain as thin wrappers so pre-0.3 code keeps
working.
"""

from __future__ import annotations

import warnings

from .config import (
    A3316_BASE,
    DMAPS_BASE,
    AtomSpec,
)
from .library import get_residue


def build_dynamic_residue(name: str, extra_bridge: int = 0) -> tuple[str, list[AtomSpec]]:
    """Look up a residue by name and return its SMILES + atom list.

    * For :class:`BetaineBase` residues, ``extra_bridge`` controls the
      number of extra CH2 groups in the head→tail linker.
    * For :class:`ResidueFragment` residues, ``extra_bridge`` must be 0.
    """
    return get_residue(name).build(extra_bridge)


# ---------------------------------------------------------------------------
# Backwards-compatible helpers (pre-0.3 API).
# ---------------------------------------------------------------------------
def build_dynamic_dmaps(extra_ch2: int) -> tuple[str, list[AtomSpec]]:
    return DMAPS_BASE.build(extra_ch2)


def build_dynamic_a3316(extra_ch2: int) -> tuple[str, list[AtomSpec]]:
    return A3316_BASE.build(extra_ch2)


def build_dynamic_betaine(name: str, extra_ch2: int) -> tuple[str, list[AtomSpec]]:
    warnings.warn(
        "build_dynamic_betaine is deprecated; use build_dynamic_residue.",
        DeprecationWarning,
        stacklevel=2,
    )
    return build_dynamic_residue(name, extra_ch2)


__all__ = [
    "build_dynamic_residue",
    "build_dynamic_dmaps",
    "build_dynamic_a3316",
    "build_dynamic_betaine",
]
