"""PDB writing + heavy/hydrogen interleaving.

Two responsibilities:

* :func:`write_pdb` — dump an RDKit molecule (labelled by :mod:`builder`) into
  a minimal PDB file, respecting the atom properties.  Emits ``CONECT``
  records so viewers render the correct bond graph regardless of geometry.
* :func:`interleave_pdb` — post-process a PDB so hydrogens follow their parent
  heavy atom (required by many MD pre-processors that expect a specific
  atom ordering per residue).  Preserves and remaps ``CONECT`` records.

Small monovalent counter-ions (Cl⁻, Br⁻, Na⁺, K⁺, …) are intentionally
**not** written here.  Add them later during solvation using Packmol or
``gmx genion`` so ionic strength matches your MD box.
"""

from __future__ import annotations

import logging
from collections.abc import Mapping

import numpy as np
from rdkit import Chem

log = logging.getLogger(__name__)

# Distance threshold (Å) below which a hydrogen is considered bonded to a
# heavy atom.  Real C-H is ~1.09 Å; O-H, N-H are slightly shorter.
_H_BOND_CUTOFF: float = 1.3


def write_pdb(mol: Chem.Mol, path: str) -> str:
    """Write a labelled molecule to ``path`` as a minimal PDB file.

    Atoms are emitted **residue by residue** (ascending residue number), each
    heavy atom immediately followed by its bonded hydrogens.  This ordering is
    required by MD pre-processors such as ``gmx pdb2gmx``, which start a new
    residue whenever the residue number changes: RDKit appends every hydrogen
    after every heavy atom, so a naive index walk would emit residue numbers
    like ``1,2,3,…,N,1,2,3,…`` and pdb2gmx would treat the repeated numbers as
    duplicate residues and fail.  The hydrogen→parent mapping uses the RDKit
    bond graph (not interatomic distance), so it is exact regardless of
    geometry — no separate distance-based interleave pass is needed.

    Emits ``CONECT`` records covering every bond in the RDKit molecule so
    downstream viewers render the correct connectivity regardless of the
    inter-atomic distances (bond guessing by distance can miss slightly
    stretched C–N⁺ bonds).
    """
    conf = mol.GetConformer()
    order = _residue_ordered_indices(mol)
    old_to_new = {old_idx: serial for serial, old_idx in enumerate(order, 1)}
    with open(path, "w") as fh:
        for serial, old_idx in enumerate(order, 1):
            atom = mol.GetAtomWithIdx(old_idx)
            pos = conf.GetAtomPosition(old_idx)
            name = (
                atom.GetProp("name")
                if atom.HasProp("name")
                else f"{atom.GetSymbol()}{serial}"
            )
            res = atom.GetProp("res") if atom.HasProp("res") else "UNK"
            resnum = int(atom.GetProp("resnum")) if atom.HasProp("resnum") else 1
            fh.write(
                f"ATOM  {serial:5d} {name:<4s} {res:3s} A{resnum:4d}    "
                f"{pos.x:8.3f}{pos.y:8.3f}{pos.z:8.3f}  1.00  0.00           "
                f"{atom.GetSymbol()}\n"
            )
        fh.write("TER\n")
        _write_conect_records(fh, mol, order, old_to_new)
        fh.write("END\n")
    return path


def _residue_ordered_indices(mol: Chem.Mol) -> list[int]:
    """Return RDKit atom indices in pdb2gmx-friendly order.

    Groups atoms by residue number (ascending) and, within each residue,
    emits every heavy atom in its original index order immediately followed
    by the hydrogens bonded to it.  Atoms without a ``resnum`` tag sort to the
    end so nothing is silently dropped.
    """
    def resnum_of(atom: Chem.Atom) -> int:
        if atom.HasProp("resnum"):
            try:
                return int(atom.GetProp("resnum"))
            except ValueError:
                pass
        return 1_000_000_000  # untagged atoms go last, preserving them

    heavy_by_res: dict[int, list[Chem.Atom]] = {}
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "H":
            continue
        heavy_by_res.setdefault(resnum_of(atom), []).append(atom)

    order: list[int] = []
    seen: set[int] = set()
    for rnum in sorted(heavy_by_res):
        for heavy in heavy_by_res[rnum]:
            order.append(heavy.GetIdx())
            seen.add(heavy.GetIdx())
            for nb in heavy.GetNeighbors():
                if nb.GetSymbol() == "H" and nb.GetIdx() not in seen:
                    order.append(nb.GetIdx())
                    seen.add(nb.GetIdx())

    # Sweep up anything not reached above (isolated H, untagged heavy atoms).
    for atom in mol.GetAtoms():
        if atom.GetIdx() not in seen:
            order.append(atom.GetIdx())
            seen.add(atom.GetIdx())
    return order


def _write_conect_records(
    fh, mol: Chem.Mol, order: list[int], old_to_new: Mapping[int, int]
) -> None:
    """Emit one CONECT line per atom listing its bonded partners.

    Serials use the post-reorder numbering (``old_to_new``) so the records
    reference the same atoms written above.  Lines are emitted in ascending
    serial order (i.e. in ``order``).
    """
    for old_idx in order:
        atom = mol.GetAtomWithIdx(old_idx)
        central = old_to_new[old_idx]
        neighbors = sorted(old_to_new[nb.GetIdx()] for nb in atom.GetNeighbors())
        if not neighbors:
            continue
        # The PDB CONECT format allows up to 4 partners per line; wrap if more.
        for chunk_start in range(0, len(neighbors), 4):
            chunk = neighbors[chunk_start:chunk_start + 4]
            parts = "".join(f"{n:>5d}" for n in chunk)
            fh.write(f"CONECT{central:>5d}{parts}\n")


def _element_of(line: str) -> str:
    """Return the element symbol from a PDB ATOM/HETATM line.

    Uses the element column (77-78) when present.  Falls back to inferring
    from the atom name only if the element column is empty — but any
    reasonable writer (including :func:`write_pdb`) populates it.
    """
    element = line[76:78].strip() if len(line) >= 78 else ""
    if element:
        return element.upper()
    # Fallback for old-style / stripped PDBs: first non-digit character.
    name = line[12:16].strip()
    for ch in name:
        if ch.isalpha():
            return ch.upper()
    return ""


def _is_hydrogen(line: str) -> bool:
    """True iff the PDB line describes a hydrogen (checked via element column)."""
    return _element_of(line) == "H"


def interleave_pdb(input_pdb: str, output_pdb: str) -> str:
    """Rewrite ``input_pdb`` so hydrogens immediately follow their heavy parent.

    A hydrogen is assigned to the closest heavy atom within
    :data:`_H_BOND_CUTOFF` Å, restricted to the same residue.  ``CONECT``
    records are also remapped to reference the new post-interleave serial
    numbers so bond graphs survive the reordering.
    """
    heavy: list[dict] = []
    hydrogens: list[dict] = []
    conect_lines: list[str] = []

    with open(input_pdb) as fh:
        for line in fh:
            if line.startswith("CONECT"):
                conect_lines.append(line)
                continue
            if not line.startswith("ATOM"):
                continue
            name = line[12:16].strip()
            res_num = line[22:26].strip()
            coords = np.array(
                [float(line[30:38]), float(line[38:46]), float(line[46:54])]
            )
            atom = {
                "line": line,
                "name": name,
                "res_num": res_num,
                "coords": coords,
                "old_serial": int(line[6:11]),
                "bonded_h": [],
            }
            if _is_hydrogen(line):
                hydrogens.append(atom)
            else:
                heavy.append(atom)

    for h in hydrogens:
        best_dist = _H_BOND_CUTOFF
        best_parent: dict | None = None
        for hv in heavy:
            if h["res_num"] != hv["res_num"]:
                continue
            d = float(np.linalg.norm(h["coords"] - hv["coords"]))
            if d < best_dist:
                best_dist = d
                best_parent = hv
        if best_parent is not None:
            best_parent["bonded_h"].append(h)
        else:
            log.warning(
                "Hydrogen %s in residue %s has no heavy parent within %.2f Å",
                h["name"], h["res_num"], _H_BOND_CUTOFF,
            )

    old_to_new: dict[int, int] = {}
    with open(output_pdb, "w") as fh:
        fh.write("TITLE     INTERLEAVED HEAVY-H PDB\n")
        idx = 1
        for hv in heavy:
            old_to_new[hv["old_serial"]] = idx
            fh.write(_reindexed(hv["line"], idx))
            idx += 1
            for h in hv["bonded_h"]:
                old_to_new[h["old_serial"]] = idx
                fh.write(_reindexed(h["line"], idx))
                idx += 1
        fh.write("TER\n")
        for line in conect_lines:
            remapped = _remap_conect(line, old_to_new)
            if remapped is not None:
                fh.write(remapped)
        fh.write("END\n")

    log.info("Interleaved PDB written to %s", output_pdb)
    return output_pdb


def _remap_conect(line: str, old_to_new: Mapping[int, int]) -> str | None:
    """Rewrite a CONECT line using the old→new serial mapping.

    Returns the new line, or None if the referenced central atom no longer
    exists (e.g. an atom that was silently dropped upstream).
    """
    try:
        central_old = int(line[6:11])
    except ValueError:
        return None
    central_new = old_to_new.get(central_old)
    if central_new is None:
        return None
    partners: list[int] = []
    for start in (11, 16, 21, 26):
        field = line[start:start + 5]
        if not field.strip():
            break
        try:
            partners.append(int(field))
        except ValueError:
            break
    new_partners = [old_to_new[p] for p in partners if p in old_to_new]
    parts = "".join(f"{n:>5d}" for n in new_partners)
    return f"CONECT{central_new:>5d}{parts}\n"


def _reindexed(line: str, new_idx: int) -> str:
    """Rewrite the serial column of a PDB ATOM/HETATM line without disturbing others."""
    return f"{line[:6]}{new_idx:>5}{line[11:]}"


__all__ = ["write_pdb", "interleave_pdb"]
