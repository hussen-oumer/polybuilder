"""User-supplied monomer library — ported from ../cpp/monomers.hpp.

Load into polybuilder with either

    polybuilder build --library examples/cpp_monomers.py --comonomer A3367 ...

or programmatically::

    import polybuilder
    polybuilder.load_user_library("examples/cpp_monomers.py")

Two styles are demonstrated:

1. Export a ``MONOMERS`` dict — the loader will call
   :func:`polybuilder.register_residue` for each entry.
2. Or call :func:`polybuilder.register_residue` yourself at module level.

Both work in the same file.

Ported monomers
---------------
* Sulfobetaines (:class:`BetaineBase`): M3295, A3367, A3361
* Neutral vinyl monomer (:class:`ResidueFragment`): NIPAM
* Styrenic dicationic (:class:`ResidueFragment`): DABCO
* Bulky neutral acrylate (:class:`ResidueFragment`): IBOA

.. warning::

    OPLS atom types for NIPAM, DABCO, and IBOA are best-guess starting
    points. Verify them against your force-field manual (or refit charges
    with a QM protocol like RESP) *before* running production MD.
"""

from polybuilder import AtomSpec, BetaineBase, Initiator, ResidueFragment
from polybuilder.config import A3316_BASE, DMAPS_BASE

# ---------------------------------------------------------------------------
# Sulfobetaines (same head+bridge+tail shape as DMAPS / A3316)
# ---------------------------------------------------------------------------

# A3367 — acrylate analogue of DMAPS (no α-methyl on the vinyl carbon).
A3367_BASE = BetaineBase(
    head=ResidueFragment(
        smiles="C(C(=O)OCC[N+](C)(C)",
        atoms=(
            AtomSpec("CC", "opls_140", 0.0000, 4),
            AtomSpec("CD", "opls_465", 0.0000, 5),
            AtomSpec("OD2", "opls_466", 0.0000, 5),
            AtomSpec("OD1", "opls_467", 0.0000, 5),
            AtomSpec("CE", "opls_490", 0.0000, 5),
            AtomSpec("C6", "opls_136", 0.0000, 3),
            AtomSpec("N1", "opls_103", 0.0000, 3),
            AtomSpec("C4", "opls_135", 0.0000, 3),
            AtomSpec("C5", "opls_135", 0.0000, 3),
        ),
    ),
    tail=DMAPS_BASE.tail,
    bridge_ff_type="opls_136",
    bridge_charge=0.0,
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
        AtomSpec("CC", "opls_140", 0.0000, 1),   # α-C from CH2=CH-
        AtomSpec("CD", "opls_235", 0.0000, 2),   # amide C
        AtomSpec("OD1", "opls_236", 0.0000, 2),  # amide =O
        AtomSpec("N1", "opls_237", 0.0000, 2),   # amide N-H (secondary)
        AtomSpec("C1", "opls_137", 0.0000, 3),   # isopropyl CH
        AtomSpec("C2", "opls_135", 0.0000, 3),   # methyl
        AtomSpec("C3", "opls_135", 0.0000, 3),   # methyl
    ),
)

# IBOA — isobornyl acrylate.  14 heavy atoms in the fragment.
# Ring system atom ordering follows the SMILES traversal RDKit produces.
IBOA_RESIDUE = ResidueFragment(
    smiles="C(C(=O)OC1CCC2(C)C(C)(C)C1C2)",
    atoms=(
        AtomSpec("CC", "opls_140", 0.0000, 1),   # α-C (acrylate)
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

# VBD — vinyl-benzyl-DABCO dication.  The residue polymerises through the
# styrenic C=C on the benzene ring; the DABCO cage stays intact as a
# non-polymerising pendant carrying two positive charges.  20 heavy atoms
# in the post-polymerisation fragment.
#
# Named "VBD" (Vinyl-Benzyl-DABCO) rather than plain "DABCO" so it's clear
# the pendant, not the DABCO cage itself, is the polymerising piece.
#
# The N⁺s carry a net +2 charge per repeat unit.  Add the balancing
# counter-ions later during solvation (Packmol / gmx genion) so ionic
# strength matches your MD box.
#
# Aromatic + quaternary-N OPLS types are placeholders — refine before MD.
VBD_RESIDUE = ResidueFragment(
    smiles="C(c1ccc(C[N+]23CC[N+](CCCC)(CC2)CC3)cc1)",
    atoms=(
        AtomSpec("CC",  "opls_141", 0.0000, 1),   # α-C (styrenic)
        AtomSpec("CG",  "opls_145", 0.0000, 2),   # aromatic ipso (Phe convention)
        AtomSpec("CD1", "opls_145", 0.0000, 2),   # aromatic (ortho)
        AtomSpec("CE1", "opls_145", 0.0000, 2),   # aromatic (meta)
        AtomSpec("CZ",  "opls_145", 0.0000, 2),   # aromatic para (subst.)
        AtomSpec("CBZ", "opls_136", 0.0000, 3),   # benzyl -CH2-N+
        AtomSpec("N1",  "opls_103", 0.0000, 3),   # DABCO N+ (benzyl side)
        AtomSpec("CN1", "opls_136", 0.0000, 3),   # bridge 1: CH2
        AtomSpec("CN2", "opls_136", 0.0000, 3),   # bridge 1: CH2
        AtomSpec("N2",  "opls_103", 0.0000, 4),   # DABCO N+ (butyl side)
        AtomSpec("CB1", "opls_136", 0.0000, 4),   # butyl α
        AtomSpec("CB2", "opls_136", 0.0000, 4),   # butyl β
        AtomSpec("CB3", "opls_136", 0.0000, 4),   # butyl γ
        AtomSpec("CB4", "opls_135", 0.0000, 4),   # butyl terminal CH3
        AtomSpec("CN3", "opls_136", 0.0000, 3),   # bridge 2: CH2
        AtomSpec("CN4", "opls_136", 0.0000, 3),   # bridge 2: CH2
        AtomSpec("CN5", "opls_136", 0.0000, 3),   # bridge 3: CH2
        AtomSpec("CN6", "opls_136", 0.0000, 3),   # bridge 3: CH2
        AtomSpec("CE2", "opls_145", 0.0000, 2),   # aromatic (meta')
        AtomSpec("CD2", "opls_145", 0.0000, 2),   # aromatic (ortho')
    ),
)


MONOMERS = {
    # sulfobetaines
    "A3367": A3367_BASE,
    "M3295": M3295_BASE,
    "A3361": A3361_BASE,
    # single-fragment monomers (extra_bridge must be 0)
    "NIPAM": NIPAM_RESIDUE,
    "IBOA":  IBOA_RESIDUE,
    "VBD":   VBD_RESIDUE,
    # legacy alias: the old name mixed up the polymerising piece (vinyl-benzyl)
    # and the pendant (DABCO cage).  Kept so pre-existing scripts still work.
    "DABCO": VBD_RESIDUE,
}


# ---------------------------------------------------------------------------
# BVD — Bis-benzyl DABCO dication.  Both N⁺ centres of the DABCO cage carry
# a -CH2-C6H4- (benzyl-phenyl) substituent; each phenyl has a vinyl group at
# the para position.  ONE vinyl polymerises into the main polymer chain
# (backbone α-C is the first heavy atom below).  The SECOND vinyl remains
# as an *unreacted* pendant -CH=CH2 on the other benzene, which is what the
# `[…]y` notation in the reference figure denotes — a reactive group that
# can later be involved in crosslinking or a second polymerisation step.
#
# 25 heavy atoms in the post-polymerisation fragment.
BVD_RESIDUE = ResidueFragment(
    smiles="C(c1ccc(C[N+]23CC[N+](Cc4ccc(C=C)cc4)(CC2)CC3)cc1)",
    atoms=(
        # ---- polymerised (main-chain) side --------------------------------
        AtomSpec("CC",   "opls_141", 0.0000, 1),   # α-C (styrenic)
        AtomSpec("CG",   "opls_145", 0.0000, 2),   # aromatic ipso (ring 1)
        AtomSpec("CD1",  "opls_145", 0.0000, 2),   # aromatic (ring 1)
        AtomSpec("CE1",  "opls_145", 0.0000, 2),   # aromatic (ring 1)
        AtomSpec("CZ",   "opls_145", 0.0000, 2),   # aromatic para (ring 1)
        AtomSpec("CBZ",  "opls_136", 0.0000, 3),   # benzyl -CH2- (ring 1 → N1)
        AtomSpec("N1",   "opls_103", 0.0000, 3),   # DABCO N+ (polymer side)
        AtomSpec("CN1",  "opls_136", 0.0000, 3),   # bridge 1
        AtomSpec("CN2",  "opls_136", 0.0000, 3),   # bridge 1
        AtomSpec("N2",   "opls_103", 0.0000, 4),   # DABCO N+ (pendant side)
        # ---- unreacted pendant side --------------------------------------
        AtomSpec("CBZ2", "opls_136", 0.0000, 4),   # benzyl -CH2- (N2 → ring 2)
        AtomSpec("CG2",  "opls_145", 0.0000, 5),   # aromatic ipso (ring 2)
        AtomSpec("CD3",  "opls_145", 0.0000, 5),   # aromatic (ring 2)
        AtomSpec("CE3",  "opls_145", 0.0000, 5),   # aromatic (ring 2)
        AtomSpec("CZ2",  "opls_145", 0.0000, 5),   # aromatic para (ring 2)
        AtomSpec("CV1",  "opls_143", 0.0000, 5),   # pendant vinyl -CH=  (sp2)
        AtomSpec("CV2",  "opls_143", 0.0000, 5),   # pendant vinyl =CH2  (sp2)
        AtomSpec("CE4",  "opls_145", 0.0000, 5),   # aromatic (ring 2)
        AtomSpec("CD4",  "opls_145", 0.0000, 5),   # aromatic (ring 2)
        # ---- DABCO cage bridges 2 & 3 (close back to N1) -----------------
        AtomSpec("CN3",  "opls_136", 0.0000, 3),   # bridge 2
        AtomSpec("CN4",  "opls_136", 0.0000, 3),   # bridge 2 (closes to N1)
        AtomSpec("CN5",  "opls_136", 0.0000, 3),   # bridge 3
        AtomSpec("CN6",  "opls_136", 0.0000, 3),   # bridge 3 (closes to N1)
        # ---- ring 1 completion -------------------------------------------
        AtomSpec("CE2",  "opls_145", 0.0000, 2),   # aromatic (ring 1)
        AtomSpec("CD2",  "opls_145", 0.0000, 2),   # aromatic (ring 1)
    ),
)

MONOMERS["BVD"] = BVD_RESIDUE


# ---------------------------------------------------------------------------
# Example user-defined initiator (in addition to the built-in KPS/AIBA/TBHP).
# ---------------------------------------------------------------------------
# APS (ammonium persulfate) shares its sulfate-radical chemistry with KPS —
# the counter-ion doesn't appear in the polymer, so the end group is the same
# -O-SO3(-).  This entry just re-registers the KPS fragments under a
# different name so recipes reading "APS" don't have to remap.
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

INITIATORS = {
    "APS": APS_INITIATOR,
}
