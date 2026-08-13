"""Convert a PubChem-style vinyl SMILES into a starting configuration entry.

The result is intentionally a *scaffold*: OPLS atom types are placeholders
(``opls_XXX``) that a chemist must review before running MD.  The purpose is
to save typing when a new monomer arrives.
"""

from __future__ import annotations

import logging

from rdkit import Chem

from ..library.config import UNKNOWN_FF_TYPE, AtomSpec
from .exceptions import InvalidSmilesError

log = logging.getLogger(__name__)

_VINYL_PATTERN = Chem.MolFromSmarts("[CH2]=[C,CH]")


def generate_config_entry(
    res_name: str, pubchem_smiles: str
) -> tuple[str | None, list[AtomSpec] | None]:
    """Return a ``(residue_smiles, atom_list)`` pair or ``(None, None)``.

    ``(None, None)`` is returned when the SMILES either fails to parse or
    doesn't contain a terminal vinyl bond (i.e. is not a polymerisable
    monomer under our simple model).
    """
    if not pubchem_smiles:
        raise InvalidSmilesError("Empty SMILES string.")

    mol = Chem.MolFromSmiles(pubchem_smiles)
    if mol is None:
        raise InvalidSmilesError(
            f"RDKit could not parse SMILES for residue '{res_name}': {pubchem_smiles!r}"
        )

    matches = mol.GetSubstructMatches(_VINYL_PATTERN)
    if not matches:
        log.warning("No terminal vinyl (CH2=C) found in %s; skipping.", res_name)
        return None, None

    term_idx, alpha_idx = matches[0]
    alpha_atom = mol.GetAtomWithIdx(alpha_idx)
    is_meth = any(
        n.GetSymbol() == "C" and n.GetDegree() == 1 and n.GetIdx() != term_idx
        for n in alpha_atom.GetNeighbors()
    )

    edit_mol = Chem.RWMol(mol)
    edit_mol.RemoveAtom(term_idx)
    res_mol = edit_mol.GetMol()

    atom_data: list[AtomSpec] = [
        AtomSpec(
            name="CC",
            ff_type="opls_139" if is_meth else "opls_140",
            charge=0.0000,
            cgnr=1,
        )
    ]

    atom_counts: dict = {}
    for atom in res_mol.GetAtoms():
        if atom.GetIdx() == 0:
            continue
        symbol = atom.GetSymbol()
        atom_counts[symbol] = atom_counts.get(symbol, 0) + 1
        atom_data.append(
            AtomSpec(
                name=f"{symbol}{atom_counts[symbol]}",
                ff_type=UNKNOWN_FF_TYPE,
                charge=0.0000,
                cgnr=1,
            )
        )

    res_smiles = Chem.MolToSmiles(res_mol)
    formatted = f"C(C)({res_smiles})" if is_meth else f"C({res_smiles})"
    return formatted, atom_data


__all__ = ["generate_config_entry"]
