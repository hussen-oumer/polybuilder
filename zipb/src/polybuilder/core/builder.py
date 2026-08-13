"""Construct a 3D polymer molecule and label its atoms with OPLS parameters.

The pipeline is:

    1. Assemble the SMILES from a repeating (MMA)_n (comonomer)_m sequence.
    2. Embed a 3D conformer and UFF-minimise.
    3. Walk the heavy atoms in the order RDKit created them and assign
       force-field name/type/charge/charge-group + residue tags.
    4. Propagate residue tags to bonded hydrogens; number hydrogens per
       residue.

The public entry point is :func:`build_polymer`.  Residues can be either
:class:`~polybuilder.config.BetaineBase` sulfobetaines (with a tunable CH2
bridge) or single-fragment :class:`~polybuilder.config.ResidueFragment`
monomers (MMA, NIPAM, DABCO, ...).
"""

from __future__ import annotations

import logging
import warnings
from collections.abc import Sequence
from dataclasses import dataclass, field

from rdkit import Chem
from rdkit.Chem import AllChem

from ..library.config import (
    BACKBONE_BRIDGE_CHARGE,
    BACKBONE_BRIDGE_TYPE,
    BACKBONE_CAP_CHARGE,
    BACKBONE_CAP_TYPE,
    HYDROGEN_CHARGE,
    HYDROGEN_FF_TYPE,
    AtomSpec,
    BetaineBase,
    Initiator,
)
from ..library import get_initiator, get_residue
from ..support.exceptions import EmbeddingError, InvalidPolymerSpecError
from ..write.pdb import write_pdb
from ..write.rtp import write_rtp

log = logging.getLogger(__name__)


_VALID_ENDS = ("both", "head", "tail", "none")


@dataclass(frozen=True)
class PolymerSpec:
    """Parameters describing one polymer to build.

    Residues
    --------
    ``comonomer`` is looked up in the polybuilder registry (see
    :func:`polybuilder.available_residues`).  ``n_mma`` MMA units and
    ``n_comonomer`` comonomer units make up one block; the block repeats
    ``repeats`` times.  If ``cap`` is true, an extra MMA is appended at
    the tail.  ``extra_bridge`` chooses how many CH2 groups sit between
    the head and tail *of a sulfobetaine comonomer*; for single-fragment
    comonomers (NIPAM, DABCO, ...) it must be 0.

    Initiator (optional)
    --------------------
    ``initiator`` is the name of a registered :class:`Initiator` (KPS,
    AIBA, TBHP, or a user-registered one).  When ``None`` (default) the
    chain ends are simple methyl caps and no initiator fragments are
    inserted anywhere.

    Placement can be controlled in three ways (highest wins):

    * ``initiator_positions``: explicit list of positions.  Each entry is
      ``"head"``, ``"tail"``, or a 0-based integer *i* meaning "insert the
      initiator between residue *i* and residue *i+1*" (mid-chain).
    * ``initiator_percent``: target loading in %.  The tool converts it
      into a number of fragments (``round(len(sequence) * pct / 100)``)
      and places them at head, tail, and evenly-spaced mid-chain positions
      in that order of preference.
    * ``ends``: legacy shorthand — ``"both"`` / ``"head"`` / ``"tail"`` /
      ``"none"``.  Ignored when ``initiator_positions`` or
      ``initiator_percent`` is given.

    Legacy aliases
    --------------
    ``betaine`` / ``n_betaine`` remain as deprecated aliases for
    ``comonomer`` / ``n_comonomer``.
    """

    comonomer: str = ""
    first_residue: str = "MMA"
    n_first: int = 0
    n_comonomer: int = 0
    repeats: int = 1
    extra_bridge: int = 0
    cap: bool = False
    initiator: str | None = None
    ends: str = "both"
    initiator_positions: list[str | int] | None = None
    initiator_percent: float | None = None
    # Deprecated aliases (kept working so older scripts don't break).
    betaine: str | None = None
    n_betaine: int | None = None
    n_mma: int | None = None

    def __post_init__(self) -> None:
        errors: list[str] = []

        # ---- resolve deprecated aliases ----------------------------------
        if self.betaine is not None:
            warnings.warn(
                "PolymerSpec(betaine=...) is deprecated; use comonomer=...",
                DeprecationWarning,
                stacklevel=3,
            )
            if not self.comonomer:
                object.__setattr__(self, "comonomer", self.betaine)
            elif self.comonomer != self.betaine:
                errors.append("comonomer and betaine given with different values")
        if self.n_betaine is not None:
            warnings.warn(
                "PolymerSpec(n_betaine=...) is deprecated; use n_comonomer=...",
                DeprecationWarning,
                stacklevel=3,
            )
            if self.n_comonomer == 0:
                object.__setattr__(self, "n_comonomer", self.n_betaine)
            elif self.n_comonomer != self.n_betaine:
                errors.append("n_comonomer and n_betaine given with different values")
        if self.n_mma is not None:
            # n_mma is the pre-generalization name for n_first (only meaningful
            # when first_residue is MMA, which was the only choice).
            if self.first_residue != "MMA":
                errors.append(
                    "n_mma is a legacy alias for n_first that only applies when "
                    "first_residue='MMA'; use n_first with a non-MMA first_residue."
                )
            elif self.n_first == 0:
                object.__setattr__(self, "n_first", self.n_mma)
            elif self.n_first != self.n_mma:
                errors.append("n_first and n_mma given with different values")

        # ---- validate ----------------------------------------------------
        if not self.first_residue:
            errors.append("first_residue name is required")
        if not self.comonomer:
            errors.append("comonomer name is required")
        if self.n_first < 0:
            errors.append("n_first must be >= 0")
        if self.n_comonomer < 0:
            errors.append("n_comonomer must be >= 0")
        if self.repeats <= 0:
            errors.append("repeats must be >= 1")
        if self.extra_bridge < 0:
            errors.append("extra_bridge must be >= 0")
        if self.n_first == 0 and self.n_comonomer == 0 and not self.cap:
            errors.append("polymer would be empty (n_first=n_comonomer=0, cap=False)")
        if self.ends not in _VALID_ENDS:
            errors.append(f"ends must be one of {_VALID_ENDS!r} (got {self.ends!r})")

        # ---- initiator placement fields --------------------------------
        if self.initiator_percent is not None:
            if not 0 <= float(self.initiator_percent) <= 100:
                errors.append(
                    f"initiator_percent must be between 0 and 100 "
                    f"(got {self.initiator_percent})"
                )
            if self.initiator is None:
                errors.append(
                    "initiator_percent requires an initiator name (set initiator=...)"
                )
        if self.initiator_positions is not None:
            for p in self.initiator_positions:
                if isinstance(p, str):
                    if p not in ("head", "tail"):
                        errors.append(
                            f"initiator_positions entry {p!r} must be 'head', "
                            f"'tail', or an integer residue index"
                        )
                elif not isinstance(p, int):
                    errors.append(
                        f"initiator_positions entry {p!r} must be a string "
                        f"('head'/'tail') or int"
                    )
            if self.initiator is None:
                errors.append(
                    "initiator_positions requires an initiator name (set initiator=...)"
                )

        if errors:
            raise InvalidPolymerSpecError("; ".join(errors))


def homopolymer(
    monomer: str,
    n: int,
    *,
    initiator: str | None = None,
    ends: str = "both",
    bridge: int = 0,
) -> PolymerSpec:
    """Build a :class:`PolymerSpec` for a homopolymer of ``n`` ``monomer`` units.

    The simple counterpart of the ``polybuilder homo`` CLI verb.  Compose it
    with :func:`build_polymer` or :func:`~polybuilder.run_pipeline`::

        from polybuilder import homopolymer, run_pipeline
        run_pipeline(homopolymer("DMAPS", 20, initiator="KPS"), output_dir="out")

    ``bridge`` adds extra CH2 groups to a sulfobetaine spacer (ignored / must be
    0 for non-sulfobetaine monomers).
    """
    return PolymerSpec(
        first_residue=monomer,
        comonomer=monomer,   # unused (n_comonomer=0) but required to be set
        n_first=n,
        n_comonomer=0,
        repeats=1,
        extra_bridge=bridge,
        initiator=initiator,
        ends=ends,
    )


def copolymer(
    monomer_a: str,
    monomer_b: str,
    n: int | None = None,
    *,
    block: bool = False,
    na: int | None = None,
    nb: int | None = None,
    repeats: int | None = None,
    initiator: str | None = None,
    ends: str = "both",
    bridge: int = 0,
) -> PolymerSpec:
    """Build a :class:`PolymerSpec` for a two-monomer copolymer.

    The simple counterpart of the ``polybuilder copo`` CLI verb.  Two ways to
    specify the layout:

    * Simple — pass ``n``.  Alternating by default (``n`` of each, ABAB...);
      with ``block=True`` it becomes one block of ``n`` A then ``n`` B.
    * Advanced — pass any of ``na`` / ``nb`` / ``repeats`` (each defaults to 1)
      for a repeating ``(na·A + nb·B)`` block.  These win over ``n``.

    ::

        copolymer("DMAPS", "DABCO", 10)                    # alternating, 10 each
        copolymer("DMAPS", "DABCO", 10, block=True)        # 10 A then 10 B
        copolymer("DMAPS", "A3367", na=3, nb=1, repeats=4) # (AAAB) x4
    """
    if any(v is not None for v in (na, nb, repeats)):
        na_r = na if na is not None else 1
        nb_r = nb if nb is not None else 1
        repeats_r = repeats if repeats is not None else 1
    elif n is not None:
        if block:
            na_r = nb_r = n
            repeats_r = 1
        else:
            na_r = nb_r = 1
            repeats_r = n
    else:
        raise InvalidPolymerSpecError(
            "copolymer needs n (simple) or na/nb/repeats (advanced)."
        )

    return PolymerSpec(
        first_residue=monomer_a,
        comonomer=monomer_b,
        n_first=na_r,
        n_comonomer=nb_r,
        repeats=repeats_r,
        extra_bridge=bridge,
        initiator=initiator,
        ends=ends,
    )


@dataclass
class PolymerResult:
    """Return value from :func:`build_polymer`."""

    pdb_path: str
    rtp_path: str
    mol: Chem.Mol
    sequence: list[str] = field(default_factory=list)


# Two-letter prefixes for building 3-letter residue codes.  Only overrides
# the default (first two letters of the residue name, uppercased).
_PREFIX_MAP: dict[str, str] = {
    "A3316": "A3",
    "A3361": "A3",
    "A3367": "A3",
    "DMAPS": "DM",
    "M3295": "M3",
    "MMA": "MC",
    "NIPAM": "NI",
    "IBOA": "IB",
    "VBD": "VB",
    "BVD": "BV",
    "DABCO": "DA",   # legacy alias for VBD
}


def _prefix_for(name: str) -> str:
    if name in _PREFIX_MAP:
        return _PREFIX_MAP[name]
    letters = "".join(c for c in name.upper() if c.isalnum())
    return (letters + "XX")[:2]


def _resolve_initiator_placements(
    spec: PolymerSpec,
    n_residues: int,
) -> tuple[bool, bool, list[int]]:
    """Return ``(cap_head, cap_tail, mid_positions)`` for the build.

    Priority: ``initiator_positions`` > ``initiator_percent`` > ``ends``.
    Mid-position integers are 0-based residue indices; a value *i* means
    "insert an initiator fragment between residue *i* and residue *i+1*".
    Only valid indices ``0 <= i <= n_residues - 2`` are kept.
    """
    if spec.initiator is None:
        return False, False, []

    # ---- explicit position list --------------------------------------
    if spec.initiator_positions is not None:
        cap_head = "head" in spec.initiator_positions
        cap_tail = "tail" in spec.initiator_positions
        max_mid = n_residues - 2
        mids: list[int] = []
        for p in spec.initiator_positions:
            if not isinstance(p, int):
                continue
            if 0 <= p <= max_mid:
                mids.append(p)
            else:
                log.warning(
                    "initiator_positions=%s: index %d is out of range "
                    "[0, %d]; skipping",
                    spec.initiator_positions, p, max_mid,
                )
        return cap_head, cap_tail, sorted(set(mids))

    # ---- percentage-driven auto placement ----------------------------
    if spec.initiator_percent is not None:
        n = round(n_residues * float(spec.initiator_percent) / 100.0)
        if n <= 0:
            return False, False, []
        if n == 1:
            return True, False, []
        if n == 2:
            return True, True, []
        # n > 2: head + tail + (n-2) mid points spread evenly in [1, n_res-1]
        mid_slots_needed = n - 2
        max_mid = n_residues - 2
        if max_mid < 0:
            return True, True, []
        # Evenly-spaced integer indices in [0, max_mid], avoiding duplicates.
        mids = []
        if mid_slots_needed >= 1 and max_mid >= 0:
            for k in range(mid_slots_needed):
                # (k + 1) / (mid_slots_needed + 1) is a fraction in (0, 1)
                idx = round((k + 1) / (mid_slots_needed + 1) * (max_mid + 1)) - 1
                idx = max(0, min(max_mid, idx))
                mids.append(idx)
        return True, True, sorted(set(mids))

    # ---- fall back to legacy `ends` ----------------------------------
    ends = spec.ends
    return (ends in ("both", "head")), (ends in ("both", "tail")), []


def _resolve_residue(name: str, extra_bridge: int) -> tuple[str, list[AtomSpec]]:
    """Look up ``name`` and return its SMILES + atoms, using ``extra_bridge``
    for betaines and 0 for single-fragment residues."""
    residue = get_residue(name)
    if isinstance(residue, BetaineBase):
        return residue.build(extra_bridge)
    return residue.build(0)


def _assemble_sequence(spec: PolymerSpec) -> list[str]:
    block = [spec.first_residue] * spec.n_first + [spec.comonomer] * spec.n_comonomer
    if not block:
        raise InvalidPolymerSpecError(
            "Cannot build polymer: n_first and n_comonomer are both zero."
        )
    sequence = block * spec.repeats
    if spec.cap:
        sequence.append(spec.first_residue)
    return sequence


def _assemble_smiles(
    sequence: Sequence[str],
    residue_smiles: dict[str, str],
    head_frag: str = "C",
    tail_frag: str = "C",
    mid_positions: Sequence[int] = (),
    mid_frag: str = "",
) -> str:
    """Assemble the polymer SMILES.

    ``head_frag`` replaces the default leading ``C`` (methyl cap) and
    ``tail_frag`` replaces the default trailing ``C``.  For each index
    *i* in ``mid_positions``, the bridge between residues *i* and *i+1*
    is replaced by ``mid_frag`` (typically an initiator's
    ``head_smiles + tail_smiles`` combined into a single spacer).
    """
    mid_set = set(mid_positions)
    parts: list[str] = [head_frag]
    for i, res in enumerate(sequence):
        parts.append(residue_smiles[res])
        if i == len(sequence) - 1:
            parts.append(tail_frag)
        elif i in mid_set:
            parts.append(mid_frag)
        else:
            parts.append("C")
    return "".join(parts)


def _embed_3d(mol: Chem.Mol) -> Chem.Mol:
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.useRandomCoords = True
    if AllChem.EmbedMolecule(mol, params) != 0:
        raise EmbeddingError("RDKit failed to embed a 3D conformer.")
    try:
        AllChem.UFFOptimizeMolecule(mol, maxIters=1000)
    except Exception as exc:  # UFF failure is non-fatal
        log.warning("UFFOptimizeMolecule failed: %s", exc)
    return mol


def _init_residue_labels(initiator: Initiator | None) -> tuple[str, str, str]:
    """Return (head_res, tail_res, mid_res) three-letter codes for initiator atoms."""
    if initiator is None:
        return "MMC", "MMC", "MMC"
    stem = "".join(c for c in initiator.name.upper() if c.isalnum())[:2]
    # Honour explicit per-initiator block-code overrides (e.g. SO4 → KPH/KPT)
    # and fall back to the two-letter stem otherwise.
    head = initiator.head_res or f"{stem}H"
    tail = initiator.tail_res or f"{stem}T"
    mid = initiator.mid_res or f"{stem}M"
    return head, tail, mid


def _label_heavy_atoms(
    mol: Chem.Mol,
    sequence: Sequence[str],
    residue_atoms: dict[str, Sequence[AtomSpec]],
    initiator: Initiator | None = None,
    cap_head: bool = False,
    cap_tail: bool = False,
    mid_positions: Sequence[int] = (),
) -> None:
    heavy = [a for a in mol.GetAtoms() if a.GetSymbol() != "H"]
    if not heavy:
        raise InvalidPolymerSpecError("Polymer produced no heavy atoms.")

    n_heavy = len(heavy)
    n_seq = len(sequence)
    head_res_label, tail_res_label, mid_res_label = _init_residue_labels(initiator)
    mid_set = set(mid_positions)

    # Use a running resnum counter so mid-chain initiators get a unique
    # residue number that falls neatly between the residues they bridge.
    resnum = 1
    h_idx = 0

    def _apply(atom, spec_atom, name, res_label, res_num):
        atom.SetProp("name", name)
        atom.SetProp("type", spec_atom.ff_type)
        atom.SetProp("charge", f"{spec_atom.charge:.4f}")
        atom.SetProp("cgnr", str(spec_atom.cgnr))
        atom.SetProp("res", res_label)
        atom.SetProp("resnum", str(res_num))
        # Stash the per-parent hydrogen charge so _label_hydrogens can give each
        # H the OPLS charge that matches its heavy-atom neighbour.
        if spec_atom.h_charge is not None:
            atom.SetProp("hq", f"{spec_atom.h_charge:.4f}")

    # ---- HEAD cap (methyl or initiator head fragment) ------------------
    if cap_head and initiator is not None:
        for spec_atom in initiator.head_atoms:
            if h_idx >= n_heavy:
                raise InvalidPolymerSpecError(
                    "Heavy-atom map exhausted labelling head initiator atoms; "
                    "initiator.head_atoms is longer than the SMILES fragment."
                )
            _apply(heavy[h_idx], spec_atom, spec_atom.name, head_res_label, resnum)
            h_idx += 1
        # Initiator head atoms form their own residue block (e.g. KPH).
        resnum += 1
    else:
        # Default methyl head cap (CAP1).  Fold it INTO the first monomer
        # residue — same residue number and residue code — so the emitted
        # PDB/RTP match the reference aminoacids.rtp layout, whose first-residue
        # blocks (DMF, MCF, ...) contain CAP1 rather than a standalone cap
        # residue.  This mirrors how the tail CAP2 is folded into the last
        # residue.  No resnum bump: CAP1 shares the first monomer's residue.
        first_label = (_prefix_for(sequence[0]) + "F")[:3]
        cap1 = heavy[h_idx]
        cap1.SetProp("name", "CAP1")
        cap1.SetProp("res", first_label)
        cap1.SetProp("type", BACKBONE_CAP_TYPE)
        cap1.SetProp("charge", f"{BACKBONE_CAP_CHARGE:.4f}")
        cap1.SetProp("cgnr", "1")
        cap1.SetProp("resnum", str(resnum))
        h_idx += 1

    # ---- Residues + inter-residue bridges (or mid-chain initiators) ---
    for i, res_name in enumerate(sequence):
        res_num_i = resnum
        if i == 0:
            suffix = "F"
        elif i == n_seq - 1:
            suffix = "L"
        else:
            suffix = "R"
        r_label = (_prefix_for(res_name) + suffix)[:3]
        cur_atoms: Sequence[AtomSpec] = residue_atoms[res_name]

        for j, spec_atom in enumerate(cur_atoms):
            if h_idx >= n_heavy:
                raise InvalidPolymerSpecError(
                    "Heavy-atom map exhausted while labelling side chains "
                    f"(residue {i + 1}/{n_seq}, atom {j + 1}/{len(cur_atoms)}). "
                    "SMILES and atom list are out of sync — check the residue's "
                    "atom list matches the SMILES heavy-atom count."
                )
            name = "CC" if j == 0 else spec_atom.name
            _apply(heavy[h_idx], spec_atom, name, r_label, res_num_i)
            h_idx += 1
        resnum += 1

        is_last = (i == n_seq - 1)
        if is_last:
            # After the last residue: either a methyl CAP2 or an initiator
            # tail fragment.
            if cap_tail and initiator is not None:
                for spec_atom in initiator.tail_atoms:
                    if h_idx >= n_heavy:
                        raise InvalidPolymerSpecError(
                            "Heavy-atom map exhausted labelling tail initiator; "
                            "initiator.tail_atoms is longer than the SMILES fragment."
                        )
                    _apply(heavy[h_idx], spec_atom, spec_atom.name,
                           tail_res_label, resnum)
                    h_idx += 1
            elif h_idx < n_heavy:
                bridge = heavy[h_idx]
                bridge.SetProp("name", "CAP2")
                bridge.SetProp("type", BACKBONE_CAP_TYPE)
                bridge.SetProp("charge", f"{BACKBONE_CAP_CHARGE:.4f}")
                bridge.SetProp("cgnr", "1")
                bridge.SetProp("res", r_label)
                bridge.SetProp("resnum", str(res_num_i))
                h_idx += 1
        elif i in mid_set and initiator is not None:
            # Mid-chain initiator: labels heavy atoms with the initiator's
            # mid_atoms in SMILES order.  All share one residue number so
            # the RTP has a single mid-initiator block per insertion.
            if initiator.mid_atoms is None:
                raise InvalidPolymerSpecError(
                    f"initiator '{initiator.name}' has no mid_atoms defined"
                )
            for spec_atom in initiator.mid_atoms:
                if h_idx >= n_heavy:
                    raise InvalidPolymerSpecError(
                        "Heavy-atom map exhausted labelling mid-chain initiator; "
                        "mid_atoms count doesn't match the mid_smiles fragment."
                    )
                _apply(heavy[h_idx], spec_atom, spec_atom.name,
                       mid_res_label, resnum)
                h_idx += 1
            resnum += 1
        else:
            # Regular inter-residue methylene bridge (one C).
            if h_idx >= n_heavy:
                raise InvalidPolymerSpecError(
                    "Heavy-atom map exhausted at inter-residue bridge; "
                    "residue atom count likely wrong."
                )
            bridge = heavy[h_idx]
            bridge.SetProp("name", "BCH2")
            bridge.SetProp("type", BACKBONE_BRIDGE_TYPE)
            bridge.SetProp("charge", f"{BACKBONE_BRIDGE_CHARGE:.4f}")
            bridge.SetProp("cgnr", "1")
            bridge.SetProp("res", r_label)
            bridge.SetProp("resnum", str(res_num_i))
            h_idx += 1

    if h_idx != n_heavy:
        log.warning(
            "Labelled %d/%d heavy atoms; %d unlabelled tail atoms remain. "
            "This usually means the residue or initiator atom list is short.",
            h_idx, n_heavy, n_heavy - h_idx,
        )


def _label_hydrogens(mol: Chem.Mol) -> None:
    # inherit residue + charge-group from bonded heavy atom, and take the
    # hydrogen's partial charge from that parent's per-parent `hq` (set from
    # AtomSpec.h_charge) — falling back to the global aliphatic-H default.
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != "H":
            continue
        neighbors = atom.GetNeighbors()
        if not neighbors:
            continue
        parent = neighbors[0]
        for key in ("res", "resnum", "cgnr"):
            if parent.HasProp(key):
                atom.SetProp(key, parent.GetProp(key))
        hq = parent.GetProp("hq") if parent.HasProp("hq") else f"{HYDROGEN_CHARGE:.4f}"
        atom.SetProp("charge", hq)

    # unique H names per residue
    res_nums = sorted(
        {
            int(a.GetProp("resnum"))
            for a in mol.GetAtoms()
            if a.HasProp("resnum")
        }
    )
    for r_num in res_nums:
        res_h = [
            a for a in mol.GetAtoms()
            if a.GetSymbol() == "H"
            and a.HasProp("resnum")
            and int(a.GetProp("resnum")) == r_num
        ]
        for idx, h in enumerate(res_h, 1):
            h.SetProp("name", f"H{idx}")
            h.SetProp("type", HYDROGEN_FF_TYPE)


def build_polymer(
    spec: PolymerSpec,
    pdb_path: str = "polymer.pdb",
    rtp_path: str = "polymer.rtp",
) -> PolymerResult:
    """Build a labelled polymer molecule and write ``.pdb`` + ``.rtp`` files."""
    sequence = _assemble_sequence(spec)

    # Reject extra_bridge > 0 when NO residue in the polymer is a
    # sulfobetaine — otherwise the parameter would silently do nothing.
    if spec.extra_bridge > 0:
        any_betaine = any(
            isinstance(get_residue(name), BetaineBase) for name in set(sequence)
        )
        if not any_betaine:
            names = ", ".join(sorted(set(sequence)))
            raise InvalidPolymerSpecError(
                f"extra_bridge={spec.extra_bridge} but none of the residues "
                f"({names}) are sulfobetaines; nothing to apply the bridge to."
            )

    # Resolve every residue type used in the sequence exactly once.
    residue_smiles: dict[str, str] = {}
    residue_atoms: dict[str, Sequence[AtomSpec]] = {}
    for name in set(sequence):
        bridge = spec.extra_bridge if isinstance(get_residue(name), BetaineBase) else 0
        smiles, atoms = _resolve_residue(name, bridge)
        residue_smiles[name] = smiles
        residue_atoms[name] = atoms

    # Resolve initiator and figure out where it goes.
    initiator: Initiator | None = None
    cap_head = cap_tail = False
    mid_positions: list[int] = []
    if spec.initiator is not None:
        initiator = get_initiator(spec.initiator)
        cap_head, cap_tail, mid_positions = _resolve_initiator_placements(
            spec, len(sequence),
        )

    head_frag = initiator.head_smiles if (initiator and cap_head) else "C"
    tail_frag = initiator.tail_smiles if (initiator and cap_tail) else "C"

    mid_frag = ""
    if initiator and mid_positions:
        if initiator.mid_smiles is None or initiator.mid_atoms is None:
            raise InvalidPolymerSpecError(
                f"initiator '{initiator.name}' has no mid_smiles/mid_atoms; "
                "mid-chain insertion is not supported for this initiator. "
                "Restrict initiator_positions to 'head'/'tail' only, or add "
                "mid_smiles + mid_atoms to the Initiator definition."
            )
        mid_frag = initiator.mid_smiles

    smiles = _assemble_smiles(
        sequence, residue_smiles, head_frag, tail_frag,
        mid_positions=mid_positions, mid_frag=mid_frag,
    )
    log.debug(
        "Assembled SMILES (%d chars): head=%s tail=%s mids=%s\n%s",
        len(smiles), cap_head, cap_tail, mid_positions, smiles,
    )

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise InvalidPolymerSpecError(
            f"RDKit could not parse the assembled SMILES: {smiles!r}"
        )

    mol = _embed_3d(mol)
    _label_heavy_atoms(
        mol, sequence, residue_atoms,
        initiator=initiator, cap_head=cap_head, cap_tail=cap_tail,
        mid_positions=mid_positions,
    )
    _label_hydrogens(mol)

    write_pdb(mol, pdb_path)
    write_rtp(mol, rtp_path)
    log.info("Wrote %s and %s (sequence length %d)", pdb_path, rtp_path, len(sequence))

    return PolymerResult(
        pdb_path=pdb_path,
        rtp_path=rtp_path,
        mol=mol,
        sequence=sequence,
    )


__all__ = [
    "PolymerSpec",
    "PolymerResult",
    "build_polymer",
    "homopolymer",
    "copolymer",
]
