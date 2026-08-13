# polybuilder

Build 3D polymer structures (PDB) and GROMACS residue topologies (RTP) from vinyl monomer SMILES, using OPLS-AA atom types.

## Features

- One-line `homo` / `copo` commands for the two everyday cases (homopolymer, copolymer).
- Built-in monomer library — no plugin file needed: sulfobetaines (DMAPS, A3316, A3361, A3367, M3295) with a tunable CH2 bridge, plus neutral/styrenic monomers (MMA, NIPAM, IBOA, DABCO/VBD, BVD).
- Radical initiator end caps and mid-chain bridges (KPS, SO4, AIBA, TBHP, APS).
- Add your own monomers and initiators via a plugin file — no package changes.
- CLI (simple verbs, full-control `build`, and interactive) and a Python API.

## Installation

```bash
pip install .          # or: pip install -e ".[dev]" for development
```

RDKit is required. If `pip install rdkit` fails, use `conda install -c conda-forge rdkit`.

## Quick start

The two everyday cases have their own one-line commands. All the common
monomers — **DMAPS, A3361, A3367, DABCO** (plus M3295, NIPAM, IBOA, VBD, BVD) —
are **built in**, so no `--library` is needed.

```bash
# Homopolymer: 20 DMAPS units
polybuilder homo DMAPS -n 20

# ...with a KPS initiator on both chain ends
polybuilder homo DMAPS -n 20 --initiator KPS

# Homopolymer of any built-in monomer
polybuilder homo A3361 -n 15
polybuilder homo DABCO -n 10
```

```bash
# Copolymer, alternating: 10 of each (ABAB...)
polybuilder copo DMAPS DABCO -n 10

# Copolymer, block: one block of 10 A then 10 B (AAA...BBB...)
polybuilder copo DMAPS DABCO --block -n 10
```

Output is named after the recipe, e.g. `homo DMAPS -n 20 --initiator KPS`
writes `DMAPS_n20_KPS.pdb` and `DMAPS_n20_KPS.rtp`. Use `-o DIR` to choose the
directory.

Shared options for `homo`/`copo`:

- `--initiator NAME` + `--ends both|head|tail|none` — radical end caps (KPS, SO4, AIBA, TBHP, APS).
- `--bridge K` — extra CH2 groups in a sulfobetaine spacer.
- `-o DIR` — output directory (default: current directory).

For uneven or repeating copolymer blocks, `copo` also accepts
`--na`/`--nb`/`--repeats` (e.g. `copo DMAPS A3367 --na 3 --nb 1 --repeats 4`).

## Advanced: the `build` command

`build` is the full-control command behind `homo`/`copo`. Everything is a flag:

```bash
# MMA + DMAPS block copolymer, KPS initiator on both ends
polybuilder build --first-residue MMA --comonomer DMAPS \
    --n-first 3 --n-comonomer 1 --repeats 4 \
    --bridge 3 --initiator KPS --ends both -o ./out
```

- `--first-residue` defaults to `MMA`. `--n-first` / `--n-comonomer` are per-block counts; `--repeats` multiplies them.
- Load extra (non-built-in) monomers with `--library FILE` (repeatable).

### Initiator placement

Three ways to place the initiator (highest priority first):

1. `--initiator-at head,tail,5,12` — explicit positions (`head`, `tail`, and/or 0-based residue indices).
2. `--initiator-percent 5` — target ~5% of residues carry an initiator, distributed head → tail → evenly mid-chain.
3. `--ends both|head|tail|none` — shorthand for cap-only cases.

Notes:
- `SO4` adds a sulfate-radical cap for sulfobetaines (DMAPS / A3361 / A3367); head/tail only.
- Only `KPS` supports mid-chain bridges by default (`KPM` residue block). AIBA/TBHP raise an error if asked to.

### Other commands

```bash
polybuilder interactive            # interactive prompt
polybuilder list-monomers          # add --library FILE to include plugins
polybuilder list-initiators
polybuilder --help
```

## Python API

The `homopolymer()` / `copolymer()` helpers are the API counterparts of the
`homo` / `copo` verbs — they return a `PolymerSpec` you hand to `build_polymer`
or `run_pipeline`:

```python
from polybuilder import homopolymer, copolymer, run_pipeline

# Homopolymer: 20 DMAPS, KPS on both ends
run_pipeline(homopolymer("DMAPS", 20, initiator="KPS"), output_dir="out")

# Alternating copolymer: 10 of each
run_pipeline(copolymer("DMAPS", "DABCO", 10), output_dir="out")

# Block copolymer, or a repeating (3·A + 1·B) block
copolymer("DMAPS", "DABCO", 10, block=True)
copolymer("DMAPS", "A3367", na=3, nb=1, repeats=4)
```

For full control, build a `PolymerSpec` directly:

```python
from polybuilder import build_polymer, PolymerSpec

spec = PolymerSpec(
    first_residue="MMA", comonomer="DMAPS",
    n_first=3, n_comonomer=1, repeats=4, extra_bridge=3,
)
build_polymer(spec, "poly.pdb", "poly.rtp")
```

## Adding your own monomers

Write a plain Python file exposing a `MONOMERS` dict; the loader picks it up automatically.

```python
# my_monomers.py
from polybuilder import AtomSpec, BetaineBase, ResidueFragment
from polybuilder.library.config import DMAPS_BASE

# Sulfobetaine: head + tunable CH2 bridge + tail
MYBET = BetaineBase(
    head=ResidueFragment(
        smiles="C(C)(C(=O)OCC[N+](C)(C)",
        atoms=(AtomSpec("CC", "opls_139", 0.0, 4), ...),  # one per heavy atom, in SMILES order
    ),
    tail=DMAPS_BASE.tail,
    bridge_ff_type="opls_136", bridge_charge=0.0, bridge_cgnr=3,
)

# Single-fragment neutral / styrenic monomer (drop the CH2=, keep the α-C fragment)
MY_NEUTRAL = ResidueFragment(
    smiles="C(C(=O)NC(C)C)",
    atoms=(AtomSpec("CC", "opls_140", 0.0, 1), ...),
)

MONOMERS = {"MYBET": MYBET, "MYNEUT": MY_NEUTRAL}
```

Each fragment needs one `AtomSpec` per heavy atom, in RDKit's SMILES traversal order. Check the order with:

```python
from rdkit import Chem
m = Chem.MolFromSmiles("C(C(=O)NC(C)C)")
for i, a in enumerate(a for a in m.GetAtoms() if a.GetSymbol() != "H"):
    print(i, a.GetSymbol())
```

A wrong count raises `InvalidPolymerSpecError: Heavy-atom map exhausted...`. For a starter block with placeholder atom types, run `polybuilder pubchem <NAME> <SMILES>`.

## Adding your own initiators

Same plugin file — expose an `INITIATORS` dict:

```python
from polybuilder import AtomSpec, Initiator

APS = Initiator(
    name="APS",
    head_smiles="S(=O)(=O)([O-])O",     # last atom bonds to the α-C
    head_atoms=(AtomSpec("S1", "opls_493", 0.0, 1), ...),
    tail_smiles="OS(=O)(=O)[O-]",        # first atom bonds to the α-C
    tail_atoms=(AtomSpec("OA", "opls_467", 0.0, 1), ...),
)

INITIATORS = {"APS": APS}
```

- Head fragment's **last** atom and tail fragment's **first** atom are the attachment points.
- `--ends head|tail|both|none` selects which fragments are used.
- Initiator atoms get residue codes `<NAME[:2]>H` / `<NAME[:2]>T` (e.g. `KPH` / `KPT`).

## Example library

`examples/cpp_monomers.py` ports the full `../cpp/monomers.hpp` library:

| Monomer | Type | Note |
|---------|------|------|
| M3295 | sulfobetaine | DMAPS with a four-carbon sulfoalkyl tail |
| A3367 | sulfobetaine | Acrylate analogue of DMAPS (no α-methyl) |
| A3361 | sulfobetaine | Alias of built-in A3316 |
| NIPAM | neutral | N-isopropylacrylamide |
| IBOA | neutral | Isobornyl acrylate |
| VBD | styrenic dication | Vinyl-benzyl DABCO; polymerises through the styrene C=C, cage stays as a +2 pendant |
| BVD | styrenic dication | Bis-benzyl DABCO; second benzyl stays as an unreacted pendant (+2) |
| DABCO | alias | Legacy name for VBD |

### Counter-ions

Polybuilder writes only the covalent polymer, so cationic residues (e.g. VBD) leave a net positive charge. Add counter-ions during solvation:


```bash
gmx pdb2gmx  -f poly.pdb -o poly.gro -p topol.top -water tip3p
gmx editconf -f poly.gro -o poly_box.gro -c -d 1.0 -bt cubic
gmx solvate  -cp poly_box.gro -cs spc216.gro -o poly_sol.gro -p topol.top
gmx grompp   -f ions.mdp -c poly_sol.gro -p topol.top -o ions.tpr
gmx genion   -s ions.tpr -o poly_ions.gro -p topol.top -pname NA -nname CL -neutral
```

⚠️ For NIPAM / IBOA / VBD / BVD the OPLS atom types are best-guess placeholders — validate them before production MD.

## Development

```bash
pip install -e ".[dev]"
pytest
ruff check src tests
mypy src
```
