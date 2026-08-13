"""Built-in extended monomer/initiator catalog.

These residues used to live in ``examples/cpp_monomers.py`` and had to be
loaded with ``--library``.  They are now shipped as built-ins so the common
targets — DMAPS, A3361, A3367, DABCO — work out of the box with no plugin
file.

The core reference monomers (MMA, DMAPS, A3316) and initiators (KPS, SO4,
AIBA, TBHP) live in :mod:`polybuilder.config`.  This module adds the rest and
is merged into the runtime registry by :mod:`polybuilder.library`.

.. warning::

    OPLS atom types / charges for the single-fragment monomers (NIPAM, IBOA,
    VBD, BVD) and the acrylamide sulfobetaine (A3361/A3316) are best-guess
    starting points.  Verify them against your force field (or refit with a
    QM/RESP protocol) *before* running production MD.
"""

from __future__ import annotations

from .config import (
    A3316_BASE,
    DMAPS_BASE,
    AtomSpec,
    BetaineBase,
    Initiator,
    Residue,
    ResidueFragment,
)

# ---------------------------------------------------------------------------
# Sulfobetaines (same head + tunable CH2 bridge + tail shape as DMAPS/A3316)
# ---------------------------------------------------------------------------

# A3367 — acrylate analogue of DMAPS (no alpha-methyl on the vinyl carbon).
A3367_BASE = BetaineBase(
    head=ResidueFragment(
        smiles="C(C(=O)OCC[N+](C)(C)",
        atoms=(
            AtomSpec("CC", "opls_135", -0.6000, 4, 0.6000),  # backbone methine >CH-
            AtomSpec("CD", "opls_465", 0.2170, 5),           # ester carbonyl C
            AtomSpec("OD2", "opls_466", -0.2918, 5),         # ester =O
            AtomSpec("OD1", "opls_467", -0.2421, 5),         # ester -O-
            AtomSpec("CE", "opls_490", 0.1934, 5, 0.0615),   # -O-CH2-
            AtomSpec("C6", "opls_136", -0.1976, 3, 0.2042),  # -CH2-N+
            AtomSpec("N1", "opls_103", 0.1573, 3),           # quaternary N+ (no H)
            AtomSpec("C4", "opls_135", -0.4018, 3, 0.2042),  # N+-CH3
            AtomSpec("C5", "opls_135", -0.4018, 3, 0.2042),  # N+-CH3
        ),
    ),
    tail=DMAPS_BASE.tail,
    bridge_ff_type="opls_136",
    bridge_charge=-0.12,
    bridge_cgnr=3,
)

# M3295 — DMAPS with a four-carbon sulfoalkyl tail.
M3295_BASE = BetaineBase(
    head=DMAPS_BASE.head,
    tail=ResidueFragment(
        smiles="CCCCS(=O)(=O)[O-])",
        atoms=(
            AtomSpec("C1", "opls_136", 0.0000, 1),
            AtomSpec("C2", "opls_136", 0.0000, 1),
            AtomSpec("C3", "opls_136", 0.0000, 1),
            AtomSpec("C4S", "opls_136", 0.0000, 1),
            AtomSpec("S1", "opls_493", 0.0000, 1),
            AtomSpec("O1", "opls_494", 0.0000, 1),
            AtomSpec("O2", "opls_494", 0.0000, 1),
            AtomSpec("O3", "opls_494", 0.0000, 1),
        ),
    ),
    bridge_ff_type="opls_136",
    bridge_charge=0.0,
    bridge_cgnr=3,
)

# A3361 — acrylamide sulfobetaine; identical to the built-in A3316 chemistry.
A3361_BASE = A3316_BASE


# ---------------------------------------------------------------------------
# Single-fragment monomers (ResidueFragment).  extra_bridge must be 0.
# ---------------------------------------------------------------------------

# NIPAM — N-isopropylacrylamide.  7 heavy atoms after CH2= is removed.
NIPAM_RESIDUE = ResidueFragment(
    smiles="C(C(=O)NC(C)C)",
    atoms=(
        AtomSpec("CC", "opls_140", 0.0000, 1),   # alpha-C from CH2=CH-
        AtomSpec("CD", "opls_235", 0.0000, 2),   # amide C
        AtomSpec("OD1", "opls_236", 0.0000, 2),  # amide =O
        AtomSpec("N1", "opls_237", 0.0000, 2),   # amide N-H (secondary)
        AtomSpec("C1", "opls_137", 0.0000, 3),   # isopropyl CH
        AtomSpec("C2", "opls_135", 0.0000, 3),   # methyl
        AtomSpec("C3", "opls_135", 0.0000, 3),   # methyl
    ),
)

# IBOA — isobornyl acrylate.  14 heavy atoms in the fragment.
IBOA_RESIDUE = ResidueFragment(
    smiles="C(C(=O)OC1CCC2(C)C(C)(C)C1C2)",
    atoms=(
        AtomSpec("CC", "opls_140", 0.0000, 1),   # alpha-C (acrylate)
        AtomSpec("CD", "opls_465", 0.0000, 2),   # ester carbonyl C
        AtomSpec("OD1", "opls_466", 0.0000, 2),  # =O
        AtomSpec("OD2", "opls_467", 0.0000, 2),  # -O-
        AtomSpec("C1", "opls_490", 0.0000, 3),   # -O-CH (bornyl)
        AtomSpec("C2", "opls_136", 0.0000, 3),   # CH2 ring
        AtomSpec("C3", "opls_136", 0.0000, 3),   # CH2 ring
        AtomSpec("C4", "opls_135", 0.0000, 3),   # ring quaternary C (has methyl)
        AtomSpec("C4M", "opls_135", 0.0000, 3),  # methyl on C4
        AtomSpec("C5", "opls_135", 0.0000, 3),   # ring quaternary C (gem-dimethyl)
        AtomSpec("C5A", "opls_135", 0.0000, 3),  # methyl
        AtomSpec("C5B", "opls_135", 0.0000, 3),  # methyl
        AtomSpec("C6", "opls_136", 0.0000, 3),   # bridgehead CH
        AtomSpec("C7", "opls_136", 0.0000, 3),   # bridge CH2
    ),
)

# VBD — vinyl-benzyl-DABCO dication.  Polymerises through the styrenic C=C;
# the DABCO cage stays intact as a non-polymerising +2 pendant.
VBD_RESIDUE = ResidueFragment(
    smiles="C(c1ccc(C[N+]23CC[N+](CCCC)(CC2)CC3)cc1)",
    atoms=(
        AtomSpec("CC",  "opls_141", 0.0000, 1),   # alpha-C (styrenic)
        AtomSpec("CG",  "opls_145", 0.0000, 2),   # aromatic ipso
        AtomSpec("CD1", "opls_145", 0.0000, 2),   # aromatic (ortho)
        AtomSpec("CE1", "opls_145", 0.0000, 2),   # aromatic (meta)
        AtomSpec("CZ",  "opls_145", 0.0000, 2),   # aromatic para (subst.)
        AtomSpec("CBZ", "opls_136", 0.0000, 3),   # benzyl -CH2-N+
        AtomSpec("N1",  "opls_103", 0.0000, 3),   # DABCO N+ (benzyl side)
        AtomSpec("CN1", "opls_136", 0.0000, 3),   # bridge 1: CH2
        AtomSpec("CN2", "opls_136", 0.0000, 3),   # bridge 1: CH2
        AtomSpec("N2",  "opls_103", 0.0000, 4),   # DABCO N+ (butyl side)
        AtomSpec("CB1", "opls_136", 0.0000, 4),   # butyl alpha
        AtomSpec("CB2", "opls_136", 0.0000, 4),   # butyl beta
        AtomSpec("CB3", "opls_136", 0.0000, 4),   # butyl gamma
        AtomSpec("CB4", "opls_135", 0.0000, 4),   # butyl terminal CH3
        AtomSpec("CN3", "opls_136", 0.0000, 3),   # bridge 2: CH2
        AtomSpec("CN4", "opls_136", 0.0000, 3),   # bridge 2: CH2
        AtomSpec("CN5", "opls_136", 0.0000, 3),   # bridge 3: CH2
        AtomSpec("CN6", "opls_136", 0.0000, 3),   # bridge 3: CH2
        AtomSpec("CE2", "opls_145", 0.0000, 2),   # aromatic (meta')
        AtomSpec("CD2", "opls_145", 0.0000, 2),   # aromatic (ortho')
    ),
)

# BVD — bis-benzyl DABCO dication with a second, unreacted styrenic pendant.
BVD_RESIDUE = ResidueFragment(
    smiles="C(c1ccc(C[N+]23CC[N+](Cc4ccc(C=C)cc4)(CC2)CC3)cc1)",
    atoms=(
        # ---- polymerised (main-chain) side --------------------------------
        AtomSpec("CC",   "opls_141", 0.0000, 1),   # alpha-C (styrenic)
        AtomSpec("CG",   "opls_145", 0.0000, 2),   # aromatic ipso (ring 1)
        AtomSpec("CD1",  "opls_145", 0.0000, 2),
        AtomSpec("CE1",  "opls_145", 0.0000, 2),
        AtomSpec("CZ",   "opls_145", 0.0000, 2),   # aromatic para (ring 1)
        AtomSpec("CBZ",  "opls_136", 0.0000, 3),   # benzyl -CH2- (ring 1 -> N1)
        AtomSpec("N1",   "opls_103", 0.0000, 3),   # DABCO N+ (polymer side)
        AtomSpec("CN1",  "opls_136", 0.0000, 3),   # bridge 1
        AtomSpec("CN2",  "opls_136", 0.0000, 3),   # bridge 1
        AtomSpec("N2",   "opls_103", 0.0000, 4),   # DABCO N+ (pendant side)
        # ---- unreacted pendant side --------------------------------------
        AtomSpec("CBZ2", "opls_136", 0.0000, 4),   # benzyl -CH2- (N2 -> ring 2)
        AtomSpec("CG2",  "opls_145", 0.0000, 5),   # aromatic ipso (ring 2)
        AtomSpec("CD3",  "opls_145", 0.0000, 5),
        AtomSpec("CE3",  "opls_145", 0.0000, 5),
        AtomSpec("CZ2",  "opls_145", 0.0000, 5),   # aromatic para (ring 2)
        AtomSpec("CV1",  "opls_143", 0.0000, 5),   # pendant vinyl -CH= (sp2)
        AtomSpec("CV2",  "opls_143", 0.0000, 5),   # pendant vinyl =CH2 (sp2)
        AtomSpec("CE4",  "opls_145", 0.0000, 5),
        AtomSpec("CD4",  "opls_145", 0.0000, 5),
        # ---- DABCO cage bridges 2 & 3 (close back to N1) -----------------
        AtomSpec("CN3",  "opls_136", 0.0000, 3),   # bridge 2
        AtomSpec("CN4",  "opls_136", 0.0000, 3),   # bridge 2 (closes to N1)
        AtomSpec("CN5",  "opls_136", 0.0000, 3),   # bridge 3
        AtomSpec("CN6",  "opls_136", 0.0000, 3),   # bridge 3 (closes to N1)
        # ---- ring 1 completion -------------------------------------------
        AtomSpec("CE2",  "opls_145", 0.0000, 2),
        AtomSpec("CD2",  "opls_145", 0.0000, 2),
    ),
)


# ---------------------------------------------------------------------------
# APS (ammonium persulfate) — same sulfate-radical chemistry as KPS.
# ---------------------------------------------------------------------------
APS_INITIATOR = Initiator(
    name="APS",
    head_smiles="S(=O)(=O)([O-])O",
    head_atoms=(
        AtomSpec("S1", "opls_493", 0.0, 1),
        AtomSpec("O1", "opls_494", 0.0, 1),
        AtomSpec("O2", "opls_494", 0.0, 1),
        AtomSpec("O3", "opls_494", 0.0, 1),
        AtomSpec("OA", "opls_467", 0.0, 1),
    ),
    tail_smiles="OS(=O)(=O)[O-]",
    tail_atoms=(
        AtomSpec("OA", "opls_467", 0.0, 1),
        AtomSpec("S1", "opls_493", 0.0, 1),
        AtomSpec("O1", "opls_494", 0.0, 1),
        AtomSpec("O2", "opls_494", 0.0, 1),
        AtomSpec("O3", "opls_494", 0.0, 1),
    ),
)


# Residues added on top of config.RESIDUE_LIBRARY.  "DABCO" is a legacy alias
# for the vinyl-benzyl-DABCO monomer VBD.
EXTRA_RESIDUES: dict[str, Residue] = {
    "A3367": A3367_BASE,
    "M3295": M3295_BASE,
    "A3361": A3361_BASE,
    "NIPAM": NIPAM_RESIDUE,
    "IBOA": IBOA_RESIDUE,
    "VBD": VBD_RESIDUE,
    "DABCO": VBD_RESIDUE,
    "BVD": BVD_RESIDUE,
}

EXTRA_INITIATORS: dict[str, Initiator] = {
    "APS": APS_INITIATOR,
}


__all__ = [
    "A3367_BASE",
    "M3295_BASE",
    "A3361_BASE",
    "NIPAM_RESIDUE",
    "IBOA_RESIDUE",
    "VBD_RESIDUE",
    "BVD_RESIDUE",
    "APS_INITIATOR",
    "EXTRA_RESIDUES",
    "EXTRA_INITIATORS",
]
