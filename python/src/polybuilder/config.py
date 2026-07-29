"""Static residue definitions.

Each residue is described as a SMILES fragment plus a list of :class:`AtomSpec`
tuples that assign OPLS-AA atom types, partial charges and charge-group numbers
in the order the heavy atoms are laid out inside the fragment.

The MMA repeat unit is fully self-contained.  DMAPS/A3316 style monomers are
built dynamically at run time: the base ``head`` and ``tail`` fragments are
joined by *n* extra CH2 groups; see :mod:`polybuilder.monomers`.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Union


@dataclass(frozen=True)
class AtomSpec:
    """Force-field description of a single heavy atom in a residue.

    Attributes
    ----------
    name:     PDB / RTP atom name (unique within a residue).
    ff_type:  OPLS-AA atom-type token (e.g. ``opls_135``).
    charge:   Partial charge in |e|; 0.0 is a valid placeholder to be refined.
    cgnr:     Charge-group number used by GROMACS to bundle atoms.
    """

    name: str
    ff_type: str
    charge: float
    cgnr: int

    def as_tuple(self) -> tuple[str, str, float, int]:
        return (self.name, self.ff_type, self.charge, self.cgnr)


@dataclass(frozen=True)
class ResidueFragment:
    """One indivisible fragment (head, tail, or a whole residue like MMA).

    Can act as a *complete* residue in its own right for single-fragment
    monomers (NIPAM, IBOA, VBD, ...): use :meth:`build` and pass
    ``extra_bridge=0``.

    Small monovalent counter-ions (Cl⁻, Br⁻, Na⁺, …) are **not** written by
    the builder — add them later during solvation with Packmol or
    ``gmx genion`` so the ionic strength matches your MD box.
    """

    smiles: str
    atoms: tuple[AtomSpec, ...]

    def build(self, extra_bridge: int = 0) -> tuple[str, list[AtomSpec]]:
        """Return ``(smiles, atoms)``.  ``extra_bridge`` must be 0."""
        if extra_bridge != 0:
            raise ValueError(
                "ResidueFragment has no bridge segment; extra_bridge must be 0 "
                f"(got {extra_bridge}). Use a BetaineBase for tunable bridges."
            )
        return self.smiles, list(self.atoms)


@dataclass(frozen=True)
class BetaineBase:
    """Two-fragment sulfobetaine: head + variable CH2 bridge + tail."""

    head: ResidueFragment
    tail: ResidueFragment
    bridge_ff_type: str = "opls_136"
    bridge_charge: float = 0.0
    bridge_cgnr: int = 3

    def with_bridge(self, extra_ch2: int) -> tuple[str, list[AtomSpec]]:
        """Return combined SMILES and atom list for ``extra_ch2`` bridge CH2 groups."""
        if extra_ch2 < 0:
            raise ValueError("extra_ch2 must be >= 0")
        smiles = self.head.smiles + ("C" * extra_ch2) + self.tail.smiles
        atoms: list[AtomSpec] = list(self.head.atoms)
        for i in range(extra_ch2):
            atoms.append(
                AtomSpec(
                    name=f"CX{i + 1}",
                    ff_type=self.bridge_ff_type,
                    charge=self.bridge_charge,
                    cgnr=self.bridge_cgnr,
                )
            )
        atoms.extend(self.tail.atoms)
        return smiles, atoms

    def build(self, extra_bridge: int = 0) -> tuple[str, list[AtomSpec]]:
        """Alias for :meth:`with_bridge` — the polymorphic entry-point."""
        return self.with_bridge(extra_bridge)


# A registered residue is either a self-contained fragment (MMA / NIPAM /
# IBOA / DABCO) or a two-part sulfobetaine base (DMAPS / A3316 / ...).
Residue = Union[ResidueFragment, BetaineBase]


def _atoms(*entries: tuple[str, str, float, int]) -> tuple[AtomSpec, ...]:
    return tuple(AtomSpec(*e) for e in entries)


# ---------------------------------------------------------------------------
# MMA (methyl methacrylate) — a self-contained residue.
# ---------------------------------------------------------------------------
MMA_RESIDUE: ResidueFragment = ResidueFragment(
    smiles="C(C)(C(=O)OC)",
    atoms=_atoms(
        ("CC", "opls_139", 0.0000, 1),
        ("CB", "opls_135", 0.0000, 2),
        ("CD", "opls_465", 0.0000, 3),
        ("OD2", "opls_466", 0.0000, 3),
        ("OD1", "opls_467", 0.0000, 3),
        ("CE", "opls_490", 0.0000, 3),
    ),
)


# ---------------------------------------------------------------------------
# DMAPS — dimethyl(methacryloyloxyethyl)ammonio propanesulfonate.
# ---------------------------------------------------------------------------
DMAPS_BASE: BetaineBase = BetaineBase(
    head=ResidueFragment(
        smiles="C(C)(C(=O)OCC[N+](C)(C)",
        atoms=_atoms(
            ("CC", "opls_139", 0.0000, 4),
            ("CB", "opls_135", 0.0000, 4),
            ("CD", "opls_465", 0.0000, 5),
            ("OD2", "opls_466", 0.0000, 5),
            ("OD1", "opls_467", 0.0000, 5),
            ("CE", "opls_490", 0.0000, 5),
            ("C6", "opls_136", 0.0000, 3),
            ("N1", "opls_103", 0.0000, 3),
            ("C4", "opls_135", 0.0000, 3),
            ("C5", "opls_135", 0.0000, 3),
        ),
    ),
    tail=ResidueFragment(
        smiles="CCCS(=O)(=O)[O-])",
        atoms=_atoms(
            ("C1", "opls_136", 0.0000, 1),
            ("C2", "opls_136", 0.0000, 1),
            ("C3", "opls_136", 0.0000, 1),
            ("S1", "opls_493", 0.0000, 1),
            ("O1", "opls_494", 0.0000, 1),
            ("O2", "opls_494", 0.0000, 1),
            ("O3", "opls_494", 0.0000, 1),
        ),
    ),
    bridge_ff_type="opls_136",
    bridge_charge=0.0,
    bridge_cgnr=3,
)


# ---------------------------------------------------------------------------
# A3316 — acrylamide-derived sulfobetaine.
# ---------------------------------------------------------------------------
A3316_BASE: BetaineBase = BetaineBase(
    head=ResidueFragment(
        smiles="C(C(=O)NCCC[N+](C)(C)",
        atoms=_atoms(
            ("CC", "opls_140", 0.0000, 4),
            ("CD", "opls_235", 0.0000, 5),
            ("OD1", "opls_236", 0.0000, 5),
            ("N2", "opls_238", 0.0000, 5),
            ("C7", "opls_136", 0.0000, 6),
            ("C8", "opls_136", 0.0000, 6),
            ("C9", "opls_136", 0.0000, 6),
            ("N1", "opls_103", 0.0000, 3),
            ("C4", "opls_135", 0.0000, 3),
            ("C5", "opls_135", 0.0000, 3),
        ),
    ),
    tail=ResidueFragment(
        smiles="CCCS(=O)(=O)[O-])",
        atoms=_atoms(
            ("C1", "opls_136", 0.0000, 1),
            ("C2", "opls_136", 0.0000, 1),
            ("C3", "opls_136", 0.0000, 1),
            ("S1", "opls_493", 0.0000, 1),
            ("O1", "opls_494", 0.0000, 1),
            ("O2", "opls_494", 0.0000, 1),
            ("O3", "opls_494", 0.0000, 1),
        ),
    ),
    bridge_ff_type="opls_136",
    bridge_charge=0.0,
    bridge_cgnr=6,
)


BETAINE_LIBRARY: dict[str, BetaineBase] = {
    "DMAPS": DMAPS_BASE,
    "A3316": A3316_BASE,
}


# Combined built-in registry (includes single-fragment residues like MMA
# alongside the sulfobetaines).
RESIDUE_LIBRARY: dict[str, Residue] = {
    "MMA": MMA_RESIDUE,
    **BETAINE_LIBRARY,
}


# ---------------------------------------------------------------------------
# INITIATOR fragments (chain-end caps derived from a radical initiator).
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Initiator:
    """A radical initiator whose fragment caps one or both chain ends
    (and, optionally, bridges two residues mid-chain).

    ``head_smiles`` is prepended to the chain SMILES; its final atom (in
    SMILES traversal order) is bonded to the α-C of the first residue.
    ``tail_smiles`` is appended to the chain SMILES; its first atom is
    bonded to the α-C of the last residue.

    ``head_atoms`` / ``tail_atoms`` describe the OPLS parameters of each
    heavy atom, in SMILES traversal order.

    Optional mid-chain insertion
    ----------------------------
    If ``mid_smiles`` and ``mid_atoms`` are set, the initiator can also
    be inserted *between* two residues.  Its **first** parsed atom bonds
    to the previous residue's α-C and its **last** parsed atom bonds to
    the next residue's α-C.  For KPS this is a persulfate bridge
    (``-O-SO2(O⁻)-O-O-SO2(O⁻)-O-``); for asymmetric or one-ended
    initiators (AIBA, TBHP) mid-insertion is not chemically meaningful
    and both fields are ``None``.
    """

    name: str
    head_smiles: str
    tail_smiles: str
    head_atoms: tuple[AtomSpec, ...]
    tail_atoms: tuple[AtomSpec, ...]
    mid_smiles: str | None = None
    mid_atoms: tuple[AtomSpec, ...] | None = None


# --- KPS (potassium persulfate → sulfate radical: -O-SO3⁻) -----------------
KPS_INITIATOR: Initiator = Initiator(
    name="KPS",
    # HEAD fragment ends with the attachment O so `head + Cα...` bonds α-C to
    # the trailing O of the sulfate ester.
    head_smiles="S(=O)(=O)([O-])O",
    head_atoms=_atoms(
        ("S1", "opls_493", 0.0, 1),
        ("O1", "opls_494", 0.0, 1),
        ("O2", "opls_494", 0.0, 1),
        ("O3", "opls_494", 0.0, 1),   # anionic O (formal -1)
        ("OA", "opls_467", 0.0, 1),   # attachment O (ester-like)
    ),
    # TAIL fragment starts with the attachment O so `...Cα + tail` bonds
    # α-C to the leading O.
    tail_smiles="OS(=O)(=O)[O-]",
    tail_atoms=_atoms(
        ("OA", "opls_467", 0.0, 1),   # attachment O
        ("S1", "opls_493", 0.0, 1),
        ("O1", "opls_494", 0.0, 1),
        ("O2", "opls_494", 0.0, 1),
        ("O3", "opls_494", 0.0, 1),
    ),
    # MID fragment = full persulfate bridge  -O-SO2(O⁻)-O-O-SO2(O⁻)-O-
    # inserted between two residues.  10 heavy atoms, net -2 charge.
    mid_smiles="OS(=O)([O-])OOS(=O)([O-])O",
    mid_atoms=_atoms(
        ("OA1", "opls_467", 0.0, 1),  # attachment O (bonds to previous α-C)
        ("S1",  "opls_493", 0.0, 1),
        ("O1",  "opls_494", 0.0, 1),
        ("O2",  "opls_494", 0.0, 1),  # anionic O
        ("OP1", "opls_467", 0.0, 1),  # peroxide bridge O
        ("OP2", "opls_467", 0.0, 1),  # peroxide bridge O
        ("S2",  "opls_493", 0.0, 1),
        ("O3",  "opls_494", 0.0, 1),
        ("O4",  "opls_494", 0.0, 1),  # anionic O
        ("OA2", "opls_467", 0.0, 1),  # attachment O (bonds to next α-C)
    ),
)

# --- AIBA (amidinium end: -C(=NH2+)(NH2)) ----------------------------------
AIBA_INITIATOR: Initiator = Initiator(
    name="AIBA",
    head_smiles="[NH2+]=C(N)",
    head_atoms=_atoms(
        ("N1", "opls_237", 0.0, 1),   # =NH2+
        ("CA", "opls_235", 0.0, 1),   # sp2 C (attachment)
        ("N2", "opls_238", 0.0, 1),   # -NH2
    ),
    tail_smiles="C(N)=[NH2+]",
    tail_atoms=_atoms(
        ("CA", "opls_235", 0.0, 1),
        ("N2", "opls_238", 0.0, 1),
        ("N1", "opls_237", 0.0, 1),
    ),
)

# --- TBHP (tert-butoxy end: -O-C(CH3)3) ------------------------------------
TBHP_INITIATOR: Initiator = Initiator(
    name="TBHP",
    head_smiles="C(C)(C)(C)O",
    head_atoms=_atoms(
        ("CT", "opls_135", 0.0, 1),   # quaternary C
        ("CM1", "opls_135", 0.0, 1),  # methyl
        ("CM2", "opls_135", 0.0, 1),  # methyl
        ("CM3", "opls_135", 0.0, 1),  # methyl
        ("OA", "opls_154", 0.0, 1),   # ether-O attachment
    ),
    tail_smiles="OC(C)(C)C",
    tail_atoms=_atoms(
        ("OA", "opls_154", 0.0, 1),
        ("CT", "opls_135", 0.0, 1),
        ("CM1", "opls_135", 0.0, 1),
        ("CM2", "opls_135", 0.0, 1),
        ("CM3", "opls_135", 0.0, 1),
    ),
)


INITIATOR_LIBRARY: dict[str, Initiator] = {
    "KPS": KPS_INITIATOR,
    "AIBA": AIBA_INITIATOR,
    "TBHP": TBHP_INITIATOR,
}


# Default charge for placeholder hydrogen atoms.
HYDROGEN_FF_TYPE: str = "opls_140"
HYDROGEN_CHARGE: float = 0.0

# Placeholder OPLS type emitted by :mod:`polybuilder.pubchem` for atoms the
# heuristic converter cannot classify.
UNKNOWN_FF_TYPE: str = "opls_XXX"

# Backbone atom names / force-field parameters used by the builder.
BACKBONE_CAP_TYPE: str = "opls_135"
BACKBONE_BRIDGE_TYPE: str = "opls_135"

# Fields that make it easy to iterate library monomers programmatically.
KNOWN_MONOMERS: tuple[str, ...] = ("MMA",) + tuple(BETAINE_LIBRARY)


__all__ = [
    "AtomSpec",
    "ResidueFragment",
    "BetaineBase",
    "Residue",
    "Initiator",
    "MMA_RESIDUE",
    "DMAPS_BASE",
    "A3316_BASE",
    "BETAINE_LIBRARY",
    "RESIDUE_LIBRARY",
    "KPS_INITIATOR",
    "AIBA_INITIATOR",
    "TBHP_INITIATOR",
    "INITIATOR_LIBRARY",
    "HYDROGEN_FF_TYPE",
    "HYDROGEN_CHARGE",
    "UNKNOWN_FF_TYPE",
    "BACKBONE_CAP_TYPE",
    "BACKBONE_BRIDGE_TYPE",
    "KNOWN_MONOMERS",
]
