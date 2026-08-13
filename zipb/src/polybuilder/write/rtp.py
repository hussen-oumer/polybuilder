"""GROMACS residue-topology (``.rtp``) writer + reorder helper.

Two responsibilities:

* :func:`write_rtp` — emit an ``.rtp`` from a labelled RDKit molecule.  Only
  one instance of each residue name is emitted, using the *first* occurrence
  of that residue in the molecule as the template.
* :func:`reorder_rtp_by_cgnr` — post-process an ``.rtp`` in place so atoms
  within each residue are sorted by charge-group number.
"""

from __future__ import annotations

import logging
import re

from rdkit import Chem

log = logging.getLogger(__name__)

_HEADER = "[ bondedtypes ]\n; junc  ang  dih  imp\n 1 1 1 1\n\n"

# Sub-section names inside a residue block; anything else is a residue name.
_SECTION_NAMES = {"atoms", "bonds", "angles", "dihedrals", "impropers", "bondedtypes"}

# Backbone link atoms.  An inter-residue bond whose *local* atom is one of these
# is an ordinary backbone connector (BCH2 → +CC / CC → -BCH2); any other local
# atom marks an initiator↔backbone junction (e.g. a KPS/SO4 sulfate ``OA`` bonded
# to the α-carbon ``CC``) whose crossing angles/dihedrals we spell out explicitly.
_BACKBONE_LINK_ATOMS = frozenset({"CC", "BCH2"})


def _offset_prefix(offset: int) -> str | None:
    """Map a residue-number offset to the GROMACS ``+``/``-`` prefix.

    Returns ``""`` for the same residue, ``"+"`` for the next and ``"-"`` for the
    previous one, or ``None`` when the atom is more than one residue away and so
    cannot be referenced from a single RTP block (pdb2gmx generates those terms
    itself from the bond graph).
    """
    if offset == 0:
        return ""
    if offset == 1:
        return "+"
    if offset == -1:
        return "-"
    return None


def _ref_name(atom: Chem.Atom, base_resnum: int) -> str | None:
    """Return ``atom``'s ``+``/``-``-prefixed name as seen from ``base_resnum``."""
    prefix = _offset_prefix(int(atom.GetProp("resnum")) - base_resnum)
    if prefix is None:
        return None
    return prefix + atom.GetProp("name")


def _write_junction_terms(fh, junctions, base_resnum: int) -> None:
    """Emit ``[ angles ]`` / ``[ dihedrals ]`` across each initiator↔backbone bond.

    ``junctions`` is a list of ``(local_atom, remote_atom)`` pairs where the bond
    ``local_atom-remote_atom`` crosses the residue boundary and ``local_atom`` is
    the initiator side (e.g. the sulfate ``OA`` bonded to the backbone ``CC``).
    For each junction we write:

    * angles  ``x-a-b`` (x on the initiator side) and ``a-b-y`` (y on the
      backbone side), and
    * proper dihedrals ``x-a-b-y`` about the junction bond,

    with every atom referenced through :func:`_ref_name` so backbone atoms carry
    the right ``+``/``-`` prefix.  Terms reaching an atom more than one residue
    away are skipped (pdb2gmx regenerates those from the bond graph).
    """
    if not junctions:
        return

    angle_lines: list[str] = []
    dih_lines: list[str] = []
    seen_a: set[tuple[str, ...]] = set()
    seen_d: set[tuple[str, str, str, str]] = set()

    for a_local, b_remote in junctions:
        a_ref = _ref_name(a_local, base_resnum)
        b_ref = _ref_name(b_remote, base_resnum)
        if a_ref is None or b_ref is None:
            continue
        a_nbrs = [n for n in a_local.GetNeighbors() if n.GetIdx() != b_remote.GetIdx()]
        b_nbrs = [n for n in b_remote.GetNeighbors() if n.GetIdx() != a_local.GetIdx()]

        # angles centred on the initiator atom: x - a_local - b_remote
        for x in a_nbrs:
            xr = _ref_name(x, base_resnum)
            if xr is None:
                continue
            key = tuple(sorted((xr, b_ref))) + (a_ref,)
            if key in seen_a:
                continue
            seen_a.add(key)
            angle_lines.append(f"  {xr.ljust(6)} {a_ref.ljust(6)} {b_ref}")

        # angles centred on the backbone atom: a_local - b_remote - y
        for y in b_nbrs:
            yr = _ref_name(y, base_resnum)
            if yr is None:
                continue
            key = tuple(sorted((a_ref, yr))) + (b_ref,)
            if key in seen_a:
                continue
            seen_a.add(key)
            angle_lines.append(f"  {a_ref.ljust(6)} {b_ref.ljust(6)} {yr}")

        # proper dihedrals about the junction bond: x - a_local - b_remote - y
        for x in a_nbrs:
            xr = _ref_name(x, base_resnum)
            if xr is None:
                continue
            for y in b_nbrs:
                yr = _ref_name(y, base_resnum)
                if yr is None:
                    continue
                key = (xr, a_ref, b_ref, yr)
                if key in seen_d or (yr, b_ref, a_ref, xr) in seen_d:
                    continue
                seen_d.add(key)
                dih_lines.append(
                    f"  {xr.ljust(6)} {a_ref.ljust(6)} {b_ref.ljust(6)} {yr}"
                )

    if angle_lines:
        fh.write("\n [ angles ]\n")
        fh.write("\n".join(angle_lines) + "\n")
    if dih_lines:
        fh.write("\n [ dihedrals ]\n")
        fh.write("\n".join(dih_lines) + "\n")


def write_rtp(mol: Chem.Mol, path: str) -> str:
    """Write an RTP file for ``mol`` (a molecule labelled by :mod:`builder`)."""
    res_names: list[str] = []
    for atom in mol.GetAtoms():
        if atom.HasProp("res"):
            name = atom.GetProp("res")
            if name not in res_names:
                res_names.append(name)

    with open(path, "w") as fh:
        fh.write(_HEADER)

        for r_name in res_names:
            atoms_in_res = [
                a for a in mol.GetAtoms()
                if a.HasProp("res") and a.GetProp("res") == r_name
            ]
            if not atoms_in_res:
                continue

            # Use the first instance of this residue as the RTP template.
            first_res_num = min(int(a.GetProp("resnum")) for a in atoms_in_res)
            template = [
                a for a in atoms_in_res
                if int(a.GetProp("resnum")) == first_res_num
            ]
            template.sort(key=lambda x: (x.GetSymbol() == "H", x.GetIdx()))

            fh.write(f"[ {r_name} ]\n [ atoms ]\n")
            for atom in template:
                fh.write(
                    f"  {atom.GetProp('name').ljust(6)} "
                    f"{atom.GetProp('type').ljust(10)} "
                    f"{atom.GetProp('charge').ljust(8)} "
                    f"{atom.GetProp('cgnr')}\n"
                )

            # ---- bonds (intra- and inter-residue) ------------------------
            # Inter-residue bonds are read straight from the molecular graph and
            # written with +/- notation, so backbone links (BCH2 → +CC,
            # CC → -BCH2) *and* initiator junctions (KPS/SO4 OA → +/-CC) are all
            # emitted — the latter was previously missing.
            fh.write("\n [ bonds ]\n")
            seen_bonds: set[tuple[str, str]] = set()
            t_idxs = {a.GetIdx() for a in template}
            junctions: list[tuple[Chem.Atom, Chem.Atom]] = []
            for atom in template:
                a_name = atom.GetProp("name")
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetIdx() in t_idxs:
                        bond = tuple(sorted((a_name, neighbor.GetProp("name"))))
                        if bond in seen_bonds:
                            continue
                        fh.write(f"  {bond[0].ljust(6)} {bond[1]}\n")
                        seen_bonds.add(bond)
                        continue
                    # Inter-residue neighbour.
                    if not neighbor.HasProp("resnum"):
                        continue
                    n_ref = _ref_name(neighbor, first_res_num)
                    if n_ref is None:
                        continue
                    key = (a_name, n_ref)
                    if key in seen_bonds:
                        continue
                    fh.write(f"  {a_name.ljust(6)} {n_ref}\n")
                    seen_bonds.add(key)
                    if a_name not in _BACKBONE_LINK_ATOMS:
                        junctions.append((atom, neighbor))

            # ---- junction angles + dihedrals (initiator ↔ backbone) ------
            # Only spelled out for the cap/junction residue (KPH/KPT/KPM, ...);
            # pdb2gmx builds the ordinary intra-residue and backbone terms.
            _write_junction_terms(fh, junctions, first_res_num)

            fh.write("\n")

    log.info("Wrote RTP to %s (%d residue templates)", path, len(res_names))
    return path


def reorder_rtp_by_cgnr(rtp_in: str, rtp_out: str) -> str:
    """Sort atoms within each residue by their charge-group number."""
    with open(rtp_in) as fh:
        lines = fh.readlines()

    header: list[str] = []
    res_data: dict[str, dict[str, list[str]]] = {}
    current_res: str | None = None

    for line in lines:
        stripped = line.strip()
        match = re.match(r"\[\s*(\w+)\s*\]", stripped)
        if match and match.group(1).lower() not in _SECTION_NAMES:
            current_res = match.group(1)
            res_data[current_res] = {"atoms": [], "other": []}
            continue

        if current_res is None:
            header.append(line.rstrip())
            continue

        if "[ atoms ]" in line:
            continue

        if stripped.startswith("["):
            res_data[current_res]["other"].append(line.rstrip())
        elif stripped and not stripped.startswith(";"):
            bucket = "atoms" if not res_data[current_res]["other"] else "other"
            res_data[current_res][bucket].append(line.rstrip())

    with open(rtp_out, "w") as fh:
        for h in header:
            fh.write(h + "\n")
        for res, content in res_data.items():
            fh.write(f"[ {res} ]\n [ atoms ]\n")
            # (cgnr, line) tuples sort correctly on the first element and
            # give mypy a homogeneous element type.
            parsed: list[tuple[int, str]] = [
                (int(line.split()[3]), line) for line in content["atoms"]
            ]
            parsed.sort(key=lambda entry: entry[0])
            for _cgnr, line in parsed:
                fh.write(line + "\n")
            for other in content["other"]:
                fh.write(other + "\n")
            fh.write("\n")

    log.info("Reordered RTP written to %s", rtp_out)
    return rtp_out


__all__ = ["write_rtp", "reorder_rtp_by_cgnr"]
