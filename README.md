# polybuilder

Build homopolymers and copolymers **from SMILES only**, using RDKit — for the
deSalSea forward-osmosis draw polymers and a downstream GROMACS workflow.

Add a monomer by typing one SMILES string. The engine finds the vinyl `C=C`,
grows the chain head-to-tail, and optionally caps the ends with an initiator
fragment (KPS, AIBA, TBHP). Python (ready to run) and C++ (contribution
scaffold) share the same algorithm.

**You only ever edit the SMILES library** (`monomers.py` / `monomers.hpp`).
Never touch the connector engine (`main.py` / `polymer_builder.*`).

## Python

```bash
pip install rdkit
cd python
python run_example.py
```

Use it in your own script:

```python
from main import PolymerBuilder as PB, to_smiles, write_structure
from monomers import MONOMERS, INITIATORS

poly = PB.cap_hydrogen(PB.homopolymer(MONOMERS["A3367"], 10))
print(to_smiles(poly))

capped = PB.add_initiator(PB.homopolymer(MONOMERS["A3367"], 10),
                          INITIATORS["KPS"], ends="both")
write_structure(capped, "A3367_KPS_10mer.pdb")
```

Add a monomer — one line in `python/monomers.py`:

```python
"MYMONO": ZwitterionicMonomer("MYMONO", "C=C(C)C(=O)OCC[N+](C)(C)CCCS(=O)(=O)[O-]"),
```

## C++

```bash
cd cpp && mkdir build && cd build
cmake -DRDKit_DIR=$CONDA_PREFIX/lib/cmake/rdkit ..
make
./polybuilder
```

Requires RDKit C++ libraries + Boost
(`conda install -c conda-forge rdkit librdkit-dev`).

## Add-ons

Optional `ext_*.py` files import the engine without editing it:
`ext_copolymer.py` (copolymer recipes) and `ext_gromacs.py` (`.pdb` + `.rtp`
export). Run `python run_extended_example.py`.

## Notes

- Handles terminal-vinyl free-radical monomers (acrylate, methacrylate,
  acrylamide, styrenic). Head-to-tail, atactic.
- Counter-ions are added at the MD box stage, so intermediate chains carry a
  net charge (expected).

