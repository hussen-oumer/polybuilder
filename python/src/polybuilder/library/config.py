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
    h_charge: Partial charge assigned to hydrogens bonded to *this* heavy atom.
              In OPLS-AA a hydrogen's charge depends on its parent (0.06 on a
              plain alkane C, 0.2042 on an ammonium N–CH2/N–CH3, 0.1625 on a
              CH2 next to a sulfonate, ...).  ``None`` falls back to the global
              :data:`HYDROGEN_CHARGE` default.  Heavy atoms with no hydrogens
              (S, O, quaternary C) leave this ``None``.
    """

    name: str
    ff_type: str
    charge: float
    cgnr: int
    h_charge: float | None = None

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
# Charges are OPLS-AA methyl-methacrylate values (from the reference MCF block
# in aminoacids.rtp).  5th tuple element = charge of the hydrogens on that heavy
# atom.  The repeat unit (with its neutral backbone CH2 bridge) is net neutral.
MMA_RESIDUE: ResidueFragment = ResidueFragment(
    smiles="C(C)(C(=O)OC)",
    atoms=_atoms(
        ("CC", "opls_139", 0.0000, 1),           # quaternary backbone C (no H)
        ("CB", "opls_135", -0.1800, 2, 0.0600),  # alpha-methyl -CH3
        ("CD", "opls_465", 0.2170, 3),           # ester carbonyl C
        ("OD2", "opls_466", -0.2918, 3),         # ester =O
        ("OD1", "opls_467", -0.2421, 3),         # ester -O-
        ("CE", "opls_490", 0.1324, 3, 0.0615),   # -O-CH3
    ),
)


# ---------------------------------------------------------------------------
# DMAPS — dimethyl(methacryloyloxyethyl)ammonio propanesulfonate.
# ---------------------------------------------------------------------------
# Charges are OPLS-AA sulfobetaine values (from the reference DMF block in
# aminoacids.rtp).  The zwitterion is net neutral: the +1 quaternary ammonium
# (N1 + its N-CH2/N-CH3 groups) balances the -1 sulfonate (S1 + O1/O2/O3).
DMAPS_BASE: BetaineBase = BetaineBase(
    head=ResidueFragment(
        smiles="C(C)(C(=O)OCC[N+](C)(C)",
        atoms=_atoms(
            ("CC", "opls_139", 0.0000, 4),           # quaternary backbone C (no H)
            ("CB", "opls_135", -0.1800, 4, 0.0600),  # alpha-methyl -CH3
            ("CD", "opls_465", 0.2170, 5),           # ester carbonyl C
            ("OD2", "opls_466", -0.2918, 5),         # ester =O
            ("OD1", "opls_467", -0.2421, 5),         # ester -O-
            ("CE", "opls_490", 0.1934, 5, 0.0615),   # -O-CH2-
            ("C6", "opls_136", -0.1976, 3, 0.2042),  # -CH2-N+
            ("N1", "opls_103", 0.1573, 3),           # quaternary N+ (no H)
            ("C4", "opls_135", -0.4018, 3, 0.2042),  # N+-CH3
            ("C5", "opls_135", -0.4018, 3, 0.2042),  # N+-CH3
        ),
    ),
    tail=ResidueFragment(
        smiles="CCCS(=O)(=O)[O-])",
        atoms=_atoms(
            ("C1", "opls_136", -0.1976, 1, 0.2042),  # N+-CH2-
            ("C2", "opls_136", -0.1200, 1, 0.0600),  # -CH2-
            ("C3", "opls_136", -0.4691, 1, 0.1625),  # -CH2-SO3
            ("S1", "opls_493", 1.4022, 1),           # sulfonate S (no H)
            ("O1", "opls_494", -0.7527, 1),          # sulfonate O
            ("O2", "opls_494", -0.7527, 1),
            ("O3", "opls_494", -0.7527, 1),
        ),
    ),
    bridge_ff_type="opls_136",
    bridge_charge=-0.12,   # tunable extra -CH2- spacer (H default 0.06)
    bridge_cgnr=3,
)


# ---------------------------------------------------------------------------
# A3316 — acrylamide-derived sulfobetaine.
# NOTE: charges are still 0.0 placeholders.  The only reference block for this
# chemistry in aminoacids.rtp (AAF/AAR/AAL) contains artifacts (e.g. a hydrogen
# with charge +0.6), so it is not a trustworthy source.  Populate h_charge/charge
# here from a clean OPLS parameterization before using A3316 in production.
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
    # Optional overrides for the 3-letter RTP residue codes of the head/tail/
    # mid cap blocks.  When ``None`` the builder derives them from the first
    # two alphanumeric characters of ``name`` (e.g. KPS → KPH/KPT/KPM).  The
    # SO4 sulfate-radical initiator sets these to KPH/KPT so its end caps share
    # the reference ``aminoacids.rtp`` block names.
    head_res: str | None = None
    tail_res: str | None = None
    mid_res: str | None = None


# --- KPS (potassium persulfate → sulfate radical: -O-SO3⁻) -----------------
KPS_INITIATOR: Initiator = Initiator(
    name="KPS",
    # HEAD fragment ends with the attachment O so `head + Cα...` bonds α-C to
    # the trailing O of the sulfate ester.
    # Charges: OPLS sulfate-ester head group, net -1 (from reference KPH block).
    head_smiles="S(=O)(=O)([O-])O",
    head_atoms=_atoms(
        ("S1", "opls_493", 0.70068, 1),
        ("O1", "opls_494", -0.42517, 1),
        ("O2", "opls_494", -0.42517, 1),
        ("O3", "opls_494", -0.42517, 1),   # anionic O (formal -1)
        ("OA", "opls_467", -0.42517, 1),   # attachment O (ester-like)
    ),
    # TAIL fragment starts with the attachment O so `...Cα + tail` bonds
    # α-C to the leading O.  Charges: net -1 (from reference KPT block).
    tail_smiles="OS(=O)(=O)[O-]",
    tail_atoms=_atoms(
        ("OA", "opls_467", -0.42517, 1),   # attachment O
        ("S1", "opls_493", 0.70068, 1),
        ("O1", "opls_494", -0.42517, 1),
        ("O2", "opls_494", -0.42517, 1),
        ("O3", "opls_494", -0.42517, 1),
    ),
    # MID fragment = full persulfate bridge  -O-SO2(O⁻)-O-O-SO2(O⁻)-O-
    # NOTE: mid-chain persulfate charges are left at 0.0 (no clean reference
    # block yet); refine before using KPS mid-chain insertion in production.
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

# --- SO4 (sulfate radical end cap: backbone-C-O-SO3⁻) -----------------------
# Same sulfate-radical chemistry as KPS (potassium persulfate → •O-SO3⁻), but
# exposed as its own selectable initiator for the DMAPS / A3361 / A3367
# sulfobetaines.  Chosen with ``--initiator SO4``; when no initiator is
# selected the chain ends stay the default CAP1 / CAP2 methyls.
#
# The head/tail caps get the reference ``aminoacids.rtp`` block names KPH / KPT
# (via head_res / tail_res) so the SO4 group cleanly stands in for the CAP1 /
# CAP2 methyls it replaces.  Head/tail only — no mid-chain persulfate bridge.
SO4_INITIATOR: Initiator = Initiator(
    name="SO4",
    # HEAD fragment ends with the attachment O so `head + Cα...` bonds the
    # α-C to the trailing ester-like O of the sulfate.
    # Charges: OPLS sulfate-ester head group, net -1 (from reference KPH block).
    head_smiles="S(=O)(=O)([O-])O",
    head_atoms=_atoms(
        ("S1", "opls_493", 0.70068, 1),
        ("O1", "opls_494", -0.42517, 1),
        ("O2", "opls_494", -0.42517, 1),
        ("O3", "opls_494", -0.42517, 1),   # anionic O (formal -1)
        ("OA", "opls_467", -0.42517, 1),   # attachment O (ester-like)
    ),
    # TAIL fragment starts with the attachment O so `...Cα + tail` bonds the
    # α-C to the leading O.  Charges: net -1 (from reference KPT block).
    tail_smiles="OS(=O)(=O)[O-]",
    tail_atoms=_atoms(
        ("OA", "opls_467", -0.42517, 1),   # attachment O
        ("S1", "opls_493", 0.70068, 1),
        ("O1", "opls_494", -0.42517, 1),
        ("O2", "opls_494", -0.42517, 1),
        ("O3", "opls_494", -0.42517, 1),
    ),
    head_res="KPH",
    tail_res="KPT",
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
    "SO4": SO4_INITIATOR,
    "AIBA": AIBA_INITIATOR,
    "TBHP": TBHP_INITIATOR,
}


# Default charge for a hydrogen whose parent heavy atom does not override it
# via ``AtomSpec.h_charge``.  0.06 |e| is the OPLS-AA value for a hydrogen on a
# plain aliphatic carbon (opls_140 on opls_135/136/139).
HYDROGEN_FF_TYPE: str = "opls_140"
HYDROGEN_CHARGE: float = 0.06

# Placeholder OPLS type emitted by :mod:`polybuilder.pubchem` for atoms the
# heuristic converter cannot classify.
UNKNOWN_FF_TYPE: str = "opls_XXX"

# Backbone atom names / force-field parameters used by the builder.  The methyl
# end cap (CAP1/CAP2) is a neutral -CH3 group (C -0.18, 3 H +0.06); the
# inter-residue methylene bridge (BCH2) is a neutral -CH2- (C -0.12, 2 H +0.06).
BACKBONE_CAP_TYPE: str = "opls_135"
BACKBONE_CAP_CHARGE: float = -0.18
BACKBONE_BRIDGE_TYPE: str = "opls_135"
BACKBONE_BRIDGE_CHARGE: float = -0.12

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
    "SO4_INITIATOR",
    "AIBA_INITIATOR",
    "TBHP_INITIATOR",
    "INITIATOR_LIBRARY",
    "HYDROGEN_FF_TYPE",
    "HYDROGEN_CHARGE",
    "UNKNOWN_FF_TYPE",
    "BACKBONE_CAP_TYPE",
    "BACKBONE_CAP_CHARGE",
    "BACKBONE_BRIDGE_TYPE",
    "BACKBONE_BRIDGE_CHARGE",
    "KNOWN_MONOMERS",
]
