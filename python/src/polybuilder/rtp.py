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

            fh.write("\n [ bonds ]\n")
            seen_bonds = set()
            t_idxs = {a.GetIdx() for a in template}
            for atom in template:
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetIdx() not in t_idxs:
                        continue
                    bond = tuple(sorted((atom.GetProp("name"), neighbor.GetProp("name"))))
                    if bond in seen_bonds:
                        continue
                    fh.write(f"  {bond[0].ljust(6)} {bond[1]}\n")
                    seen_bonds.add(bond)

            # Standard backbone connections used by pdb2gmx.
            if not r_name.endswith("F") and r_name != "MMC":
                fh.write("  CC      -BCH2\n")
            if not r_name.endswith("L") and r_name != "MMC":
                fh.write("  BCH2    +CC\n")
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
