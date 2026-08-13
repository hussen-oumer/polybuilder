# polybuilder


![desalsea](https://desalsea.eu/) · ![m-era.net](https://www.m-era.net/joint-calls/joint-call-2024) · ![Dipc](https://scc.dipc.org/docs/systems/hyperion)

This project is developed for `[desalsea](https://desalsea.eu/)` with funding from the `[M-ERA.NET Joint Call 2024](https://www.m-era.net/joint-calls/joint-call-2024)`, and simulated on the `[DIPC Hyperion](https://scc.dipc.org/docs/systems/hyperion)` HPC system.

Build 3D polymer structures (PDB) and GROMACS residue topologies (RTP) from
 monomer SMILES, using OPLS-AA atom types.

One-line `homo` / `copo` commands cover the everyday cases, common monomers
(DMAPS, A3361, A3367, DABCO, …) are built in, and you can add your own via a
plugin file — no package changes.

```bash
cd zipb
pip install .

polybuilder homo DMAPS -n 20 --initiator KPS       # homopolymer
polybuilder copo DMAPS DABCO -n 10                 # alternating copolymer
```

## Layout

- **`zipb/`** — the maintained package (CLI + Python API). See
  [`zipb/README.md`](zipb/README.md) for full usage.
- **`cpp/`** — C++ port of the SMILES connector engine (contribution scaffold).
