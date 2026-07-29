# polybuilder

Build **homopolymers and copolymers from SMILES only**, using RDKit — designed
for the deSalSea forward-osmosis draw polymers (zwitterionic sulfobetaines,
dicationic DABCO) and meant to feed a **GROMACS** workflow later.

You add a monomer by typing **one SMILES string**. The engine finds the
polymerisable vinyl bond, grows the chain head-to-tail, and (optionally) caps
the chain ends with an **initiator-derived fragment** (KPS sulfate, AIBA
amidinium, TBHP). Both a **Python** version (run it now) and a **C++** version
(expand it, then contribute to the GROMACS community) are included and share
the exact same algorithm.

---

## Why it is split this way

> *"main is the connector and nobody touches it; the SMILES section is user input."*

| You edit | You never touch |
|---|---|
| `python/monomers.py` / `cpp/monomers.hpp` — the SMILES library | `python/main.py` / `cpp/polymer_builder.*` — the connector engine |

Adding a new monomer = adding one line to the *monomers* file. The connector
never hard-codes any chemistry, so it keeps working for whatever you add.

---

## The algorithm (how one SMILES becomes a polymer)

A vinyl monomer polymerises through its terminal `C=C`. The engine:

1. **Finds** the terminal vinyl with SMARTS `[CH2]=[CX3]`.
2. **Opens** that double bond to a single bond and tags the two backbone
   carbons with dummy attachment points:
   - `head` = the `=CH2` carbon  (atom map `:1`)
   - `tail` = the substituted carbon (atom map `:2`)

   Everything else — the **cationic centre** `[N+]`, the **CH2 bridge**, the
   **anionic sulfonate** `S(=O)(=O)[O-]` — is left untouched. That is your
   "separate the SMILES into cationic / anionic / bridge" idea: the side chain
   is carried along as-is; only the backbone bond is edited.
3. **Joins** `tail(:2)` of unit *i* to `head(:1)` of unit *i+1*, giving the
   real vinyl backbone `-CH2-CR-CH2-CR-`.
4. **Caps** the two open ends with either H or an **initiator fragment**
   (e.g. KPS `[*]OS(=O)(=O)[O-]`, the `O-SO3-` end you described).

Because step 2 only ever edits `C=C`, **any** side chain works.

```
   C=C(C)C(=O)OCC[N+](C)(C)CCCS(=O)(=O)[O-]        (DMAPS monomer)
        │ open C=C, tag backbone
        ▼
   [*:1]-CH2-C(:2*)(CH3)-C(=O)O-CH2CH2-N+(CH3)2-CH2CH2CH2-SO3-
        │ join tail->head, ×n, cap ends
        ▼
   CC(CC(C)(...)C(=O)O...SO3-)...   (homopolymer)
```

---

## Python (use it right now)

### Install
```bash
pip install rdkit
```

### Run the demo
```bash
cd python
python run_example.py
```

### Use it in your own script
```python
from main import PolymerBuilder as PB
from main import to_smiles, write_structure
from monomers import MONOMERS, INITIATORS

# homopolymer of 10 A3367 units
poly = PB.cap_hydrogen(PB.homopolymer(MONOMERS["A3367"], 10))
print(to_smiles(poly))

# copolymers
PB.alternating(MONOMERS["DMAPS"], MONOMERS["MMA"], n_pairs=5)
PB.block(MONOMERS["DMAPS"], 8, MONOMERS["MMA"], 2)
PB.random_copolymer(MONOMERS["DABCO"], MONOMERS["NIPAM"], total=20, frac_a=0.5, seed=1)

# add the initiator end group (KPS sulfate on both ends)
chain = PB.homopolymer(MONOMERS["A3367"], 10)
capped = PB.add_initiator(chain, INITIATORS["KPS"], ends="both")

# 3D structure for the force-field step
write_structure(capped, "A3367_KPS_10mer.pdb")   # -> LigParGen / CGenFF / antechamber
```

### Add a new monomer (the only file you edit)
Open `python/monomers.py` and add a line:
```python
"MYMONO": ZwitterionicMonomer("MYMONO", "C=C(C)C(=O)OCC[N+](C)(C)CCCS(=O)(=O)[O-]"),
```
Choose the class for what it is: `ZwitterionicMonomer`, `DicationicMonomer`,
or `NeutralMonomer`. (They inherit everything; the class only labels the
chemistry for `analyze_fragments()`.)

---

## C++ (expand it, then upstream to GROMACS)

Same classes, same algorithm — the scaffold you grow.

```bash
cd cpp
mkdir build && cd build
cmake -DRDKit_DIR=$CONDA_PREFIX/lib/cmake/rdkit ..
make
./polybuilder
```
Requires the **RDKit C++ libraries** and Boost (easiest via
`conda install -c conda-forge rdkit librdkit-dev`). Edit `cpp/monomers.hpp` to
add monomers; leave `cpp/polymer_builder.*` alone.

---

## The initiator model

You said: *"if I want to add initiator I will connect it based on the input
percentage in the experimental work."* Two separate things happen:

1. **The end group** — the initiator fragment (KPS → `-O-SO3-`, AIBA →
   amidinium, TBHP → `-OtBu`) is bonded onto the chain ends by
   `add_initiator(...)`. This is the chemistry that changes solubility and
   cloud point.
2. **The loading (wt %)** — in real free-radical polymerisation the initiator
   percentage sets the **number of chains**, i.e. the molar mass, *not* the
   chemistry of one chain. `PolymerBuilder.estimate_dp(monomer_MW,
   initiator_MW, wt_percent)` turns your lab loading into a rough
   number-average degree of polymerisation, which you pass as `n`.

> `estimate_dp` is a **heuristic** (tune `radical_efficiency` / `conversion`
> to your GPC). Anchor it to the deck's measured `Mn` before trusting it.

---

## From here to GROMACS

`write_structure(...)` gives a 3D `.pdb`/`.mol`. Continue the pipeline:

1. **Parameterise** each unique molecule → LigParGen (OPLS), CGenFF (CHARMM),
   or `antechamber` (GAFF). Apply ECC charge scaling (~0.75–0.8) for the
   charged backbones.
2. **Assemble** the box → `gmx insert-molecules`, `gmx solvate`, `gmx genion`
   (add the DABCO counter-ions / KPS `K+` here).
3. **Equilibrate & run** → EM / NVT / NPT / production.
4. **Analyse** → phase behaviour (`gmx gyrate`, `clustsize`, RDF, chain–chain
   PMF) and osmotic pressure (virtual semipermeable wall).

(See the companion *deSalSea GROMACS workflow* document for the full protocol.)

---

## Add-ons (build on top, touch nothing)

The engine is meant to be **extended without editing it**. Extra features live
in their own `ext_*.py` files that only *import* the engine. Delete an add-on
and the core still runs. Two are included:

### 1. Copolymer recipes — `ext_copolymer.py`
Describe a copolymer instead of hand-writing a sequence. Any comonomer (e.g.
**MMA**) can be mixed into any base, with full control of **how many units** (`n`)
and the arrangement:

```python
from ext_copolymer import CopolymerRecipe

# 30% MMA copolymerised into DMAPS, 20 units total, random placement
recipe = CopolymerRecipe(base="DMAPS", comonomer="MMA",
                         n=20, frac_comonomer=0.3, mode="random")
seq = recipe.sequence()          # -> list of Monomer objects
```
`mode` = `random` · `alternating` · `block` · `gradient`.
`n` is the repeat count (degree of polymerisation); `frac_comonomer` the target
mole fraction. `comonomer=None` gives a homopolymer of length `n`.

### 2. GROMACS export — `ext_gromacs.py`
Writes a matching **`.pdb` + `.rtp`** with three residue types, exactly as
requested — **first residue, repeat residue, last residue**:

```python
from ext_gromacs import export_polymer
from monomers import INITIATORS

rp = export_polymer(recipe.sequence(), "DMAPS_co_MMA",
                    initiator=INITIATORS["KPS"], ends="both")
# -> DMAPS_co_MMA.pdb  and  DMAPS_co_MMA.rtp
```
The `.rtp` gives one block per distinct residue, e.g. for a DMAPS-co-MMA chain:

| resname | role | notes |
|---|---|---|
| `DMH` | first  | DMAPS head + initiator cap |
| `DMM` | repeat | DMAPS middle unit |
| `MMM` | repeat | MMA middle unit |
| `DMT` | last   | DMAPS tail + initiator cap |

Backbone residues are linked with the standard GROMACS `CB  +CA` notation
(head carbon = `CA`, tail carbon = `CB` in every residue). Atom **names**,
per-residue **bonds** and the inter-residue link are real and consistent
between the `.pdb` and `.rtp`.

> **Placeholders you must fill:** atom **charges** are Gasteiger estimates and
> atom **types** are element-based placeholders (`CT/HC/OS/NT/SO`). Replace the
> type column with your force field's atom types (GAFF/CGenFF/OPLS) before
> `pdb2gmx`. Edit the `_TYPE_GUESS` map in `ext_gromacs.py` to change the guess.

Run the add-on demo:
```bash
cd python && python run_extended_example.py
```

### Writing your own add-on
Create a new `ext_myfeature.py`, `from main import ...`, and build. Never edit
`main.py` / `monomers.py`. That keeps upgrades painless and your features
composable.

---

## Project layout
```
polybuilder/
├── README.md
├── python/
│   ├── main.py                  # connector engine — do not edit
│   ├── monomers.py              # SMILES library — EDIT THIS
│   ├── run_example.py           # core demo
│   ├── ext_copolymer.py         # ADD-ON: copolymer recipes
│   ├── ext_gromacs.py           # ADD-ON: .pdb + .rtp (first/repeat/last)
│   └── run_extended_example.py  # add-on demo
└── cpp/
    ├── polymer_builder.hpp / .cpp   # connector engine — do not edit
    ├── monomers.hpp                 # SMILES library — EDIT THIS
    ├── main.cpp                     # runnable demo
    └── CMakeLists.txt
```

## Class overview (both languages)
```
Monomer (base: reads SMILES, finds vinyl, makes repeat unit)
 ├── ZwitterionicMonomer   (DMAPS, M3295, A3367, A3361)
 ├── DicationicMonomer     (DABCO)
 └── NeutralMonomer        (MMA, NIPAM, IBOA)

Initiator            (cap fragment [*]... + experimental wt%)
PolymerBuilder       (homopolymer / copolymer / alternating / block /
                      random_copolymer / add_initiator / estimate_dp)
```

## Limitations / notes
- Handles **terminal-vinyl** free-radical monomers (acrylate, methacrylate,
  acrylamide, styrenic). Ring-opening / condensation monomers are out of scope.
- Head-to-tail addition, atactic (no stereochemistry set) — fine for building
  MD topologies; add tacticity later if needed.
- DABCO counter-ions (Cl⁻/Br⁻) are added at the MD box stage, not in the
  backbone SMILES, so intermediate chains carry a net positive charge (expected).
- The Python engine is tested end-to-end; the C++ engine mirrors it line-for-line
  and is the intended contribution scaffold.

## License
MIT (suggested) — add a `LICENSE` file before publishing.
