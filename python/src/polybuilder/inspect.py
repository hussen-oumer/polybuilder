"""Read a PDB and dump per-residue atoms + connectivity.

Useful for verifying a build without relying on a 3D viewer's ``CONECT``
support.  Parses the file the same way any downstream tool would: it
groups atoms by ``(chain, resnum, resname)`` and consults the file's
``CONECT`` records — no chemistry inference.

Public entry points:

* :func:`inspect_pdb`   — parse into an :class:`InspectReport` (structured).
* :func:`print_report`  — pretty-print a report to stdout.
* :func:`inspect`       — convenience wrapper (parse + print).
"""

from __future__ import annotations

import sys
from collections import defaultdict
from dataclasses import dataclass, field


@dataclass(frozen=True)
class Atom:
    serial: int
    name: str
    element: str
    resname: str
    chain: str
    resnum: int
    is_hetatm: bool


@dataclass
class Residue:
    resname: str
    chain: str
    resnum: int
    atoms: list[Atom] = field(default_factory=list)

    @property
    def heavy(self) -> list[Atom]:
        return [a for a in self.atoms if a.element != "H"]

    @property
    def hydrogens(self) -> list[Atom]:
        return [a for a in self.atoms if a.element == "H"]

    @property
    def key(self) -> tuple[str, int, str]:
        return (self.chain, self.resnum, self.resname)


@dataclass
class InspectReport:
    path: str
    atoms: dict[int, Atom]
    bonds: dict[int, set[int]]
    residues: list[Residue]

    @property
    def n_bonds(self) -> int:
        # Each undirected bond appears twice in the symmetric bonds dict.
        return sum(len(nbrs) for nbrs in self.bonds.values()) // 2


def _parse_pdb(path: str) -> tuple[dict[int, Atom], dict[int, set[int]]]:
    atoms: dict[int, Atom] = {}
    bonds: dict[int, set[int]] = defaultdict(set)
    with open(path) as fh:
        for line in fh:
            if line.startswith(("ATOM  ", "HETATM")):
                try:
                    serial = int(line[6:11])
                    name = line[12:16].strip()
                    resname = line[17:20].strip()
                    chain = line[21:22].strip() or " "
                    resnum = int(line[22:26])
                except ValueError:
                    continue
                element = line[76:78].strip() if len(line) >= 78 else ""
                if not element:
                    element = name[0]
                atoms[serial] = Atom(
                    serial=serial, name=name, element=element,
                    resname=resname, chain=chain, resnum=resnum,
                    is_hetatm=line.startswith("HETATM"),
                )
            elif line.startswith("CONECT"):
                try:
                    central = int(line[6:11])
                except ValueError:
                    continue
                for start in (11, 16, 21, 26, 31):
                    field = line[start:start + 5] if len(line) > start else ""
                    if not field.strip():
                        break
                    try:
                        nbr = int(field)
                    except ValueError:
                        break
                    bonds[central].add(nbr)
                    bonds[nbr].add(central)
    return atoms, dict(bonds)


def inspect_pdb(path: str) -> InspectReport:
    """Parse ``path`` and return a structured :class:`InspectReport`."""
    atoms, bonds = _parse_pdb(path)
    groups: dict[tuple[str, int, str], Residue] = {}
    for atom in atoms.values():
        key = (atom.chain, atom.resnum, atom.resname)
        if key not in groups:
            groups[key] = Residue(resname=atom.resname, chain=atom.chain, resnum=atom.resnum)
        groups[key].atoms.append(atom)
    residues = sorted(groups.values(), key=lambda r: (r.chain, r.resnum))
    return InspectReport(path=path, atoms=atoms, bonds=bonds, residues=residues)


def print_report(
    report: InspectReport,
    *,
    verbose: bool = False,
    residue_filter: str | None = None,
    file=sys.stdout,
) -> None:
    """Pretty-print an :class:`InspectReport`.

    * Default: one section per residue listing its heavy atoms + any bonds
      that cross into another residue (which is what you usually want to
      check — is the pendant actually attached?).
    * ``verbose=True``: also list every neighbour of every heavy atom.
    * ``residue_filter``: only include residues whose ``resname`` matches.
    """
    print(f"File: {report.path}", file=file)
    print(f"  atoms   : {len(report.atoms)}", file=file)
    print(f"  bonds   : {report.n_bonds}   (from CONECT records)", file=file)
    print(f"  residues: {len(report.residues)}", file=file)
    print(file=file)

    kept = [r for r in report.residues if residue_filter is None or r.resname == residue_filter]
    if residue_filter and not kept:
        print(f"No residues named {residue_filter!r} found.", file=file)
        return

    for res in kept:
        heavy = res.heavy
        hn = len(res.hydrogens)
        print(
            f"Residue {res.resname}  chain={res.chain}  resnum={res.resnum}  "
            f"({len(heavy)} heavy + {hn} H)",
            file=file,
        )
        print(f"  Heavy atoms: {', '.join(a.name for a in heavy)}", file=file)

        # Inter-residue bonds: quick answer to "is X actually attached to Y?".
        inter_lines: list[str] = []
        for atom in heavy:
            crossings: list[str] = []
            for nbr_serial in sorted(report.bonds.get(atom.serial, set())):
                nbr = report.atoms.get(nbr_serial)
                if nbr is None:
                    continue
                if (nbr.chain, nbr.resnum) != (res.chain, res.resnum):
                    crossings.append(f"{nbr.name}({nbr.resname}{nbr.resnum}/{nbr.chain})")
            if crossings:
                inter_lines.append(f"    {atom.name:<5s} → {', '.join(crossings)}")

        if inter_lines:
            print("  Inter-residue bonds:", file=file)
            for line in inter_lines:
                print(line, file=file)
        else:
            print("  Inter-residue bonds: none", file=file)

        if verbose:
            print("  Intra-residue bonds:", file=file)
            for atom in heavy:
                intra: list[str] = []
                for nbr_serial in sorted(report.bonds.get(atom.serial, set())):
                    nbr = report.atoms.get(nbr_serial)
                    if nbr is None:
                        continue
                    if (nbr.chain, nbr.resnum) == (res.chain, res.resnum):
                        intra.append(nbr.name)
                if intra:
                    print(f"    {atom.name:<5s} → {', '.join(intra)}", file=file)
        print(file=file)


def inspect(
    path: str,
    *,
    verbose: bool = False,
    residue_filter: str | None = None,
    file=sys.stdout,
) -> InspectReport:
    """Convenience: parse ``path`` and print the report immediately."""
    report = inspect_pdb(path)
    print_report(report, verbose=verbose, residue_filter=residue_filter, file=file)
    return report


__all__ = [
    "Atom",
    "Residue",
    "InspectReport",
    "inspect_pdb",
    "print_report",
    "inspect",
]
