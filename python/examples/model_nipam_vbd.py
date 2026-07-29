"""NIPAM + vinyl-benzyl-DABCO (VBD) copolymer with KPS sulfate end groups.

Run:  python examples/model_nipam_vbd.py

Chemistry
---------
* NIPAM polymerises through its acrylamide C=C  → -CH2-CH(CONHiPr)- backbone
* VBD    polymerises through its styrenic  C=C  → -CH2-CH(Ar-CH2-DABCO²⁺-)-
  The DABCO cage stays intact as a pendant carrying two positive charges.
* Two counter-ions per VBD unit (one Cl⁻ + one Br⁻) are placed in the PDB
  near each of the two N⁺ centers.  Chain ID for the ions is "I".
* KPS sulfate ends on both sides of the chain.

Output: out_nipam_vbd/poly.pdb  and  out_nipam_vbd/rearranged_polymer.rtp
"""

from pathlib import Path

import polybuilder

# 1. Register the extra monomers (NIPAM, VBD) + user initiators (APS).
polybuilder.load_user_library(Path(__file__).with_name("cpp_monomers.py"))

# 2. Recipe — tweak these to sweep composition.
spec = polybuilder.PolymerSpec(
    first_residue="NIPAM", n_first=3,
    comonomer="DABCO",       n_comonomer=1,
    repeats=4,
    cap=True,
    initiator="KPS", ends="both",
)

# 3. Build + write everything into ./out_nipam_vbd/.
result = polybuilder.run_pipeline(
    spec,
    output_dir="out_nipam_vbd",
    run_editconf=False,
)

n_residues = (spec.n_first + spec.n_comonomer) * spec.repeats + (1 if spec.cap else 0)
print(f"PDB: {result.pdb_path}")
print(f"RTP: {result.rtp_path}")
print(f"Sequence: ({spec.first_residue} × {spec.n_first} + "
      f"{spec.comonomer} × {spec.n_comonomer}) × {spec.repeats}"
      + (" + cap" if spec.cap else "")
      + f" — {n_residues} residues total")
