"""Minimal build script — edit the values below, then:  python model.py

Produces a PDB + RTP in ./out/ (created if needed). No CLI arguments needed.
For a fully-flagged version see the CLI: `polybuilder build --help`.
"""

from pathlib import Path

import polybuilder

# ---------------------------------------------------------------------------
# 1. Load extra monomers/initiators from a plugin file (optional).
#    Skip this line if you only need built-ins (MMA/DMAPS/A3316 + KPS/AIBA/TBHP).
# ---------------------------------------------------------------------------
polybuilder.load_user_library(Path(__file__).with_name("cpp_monomers.py"))

# ---------------------------------------------------------------------------
# 2. Describe the polymer you want to build.
#    - first_residue + comonomer :  the two monomers in each block
#    - n_first + n_comonomer     :  how many of each per block
#    - repeats                   :  how many times the block repeats
#    - cap                       :  add one extra first_residue at the tail
#    - initiator + ends          :  optional radical initiator on head/tail/both
# ---------------------------------------------------------------------------
spec = polybuilder.PolymerSpec(
    first_residue="DMAPS", n_first=1,
    comonomer="DMAPS",    n_comonomer=1,
    repeats=5,
    cap=True,
    initiator="KPS",
    ends="both",
)

# ---------------------------------------------------------------------------
# 3. Run the full pipeline (build + interleave + reorder + optional gmx step).
# ---------------------------------------------------------------------------
result = polybuilder.run_pipeline(
    spec,
    output_dir="out",
    run_editconf=False,   # set True if `gmx` is on your PATH
)

print(f"PDB: {result.pdb_path}")
print(f"RTP: {result.rtp_path}")
