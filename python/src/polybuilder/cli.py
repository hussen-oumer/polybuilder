"""Command-line entry point for polybuilder.

Two subcommands:

* ``polybuilder build`` — non-interactive, argparse-driven.
* ``polybuilder interactive`` — reproduces the original ``input()`` prompt.

Both delegate to :func:`polybuilder.pipeline.run_pipeline`.
"""

from __future__ import annotations

import argparse
import logging
import sys
from collections.abc import Sequence

from .builder import PolymerSpec
from .exceptions import PolybuilderError
from .inspect import inspect as inspect_pdb_cmd
from .library import available_initiators, available_residues, load_user_library
from .pipeline import run_pipeline
from .pubchem import generate_config_entry

log = logging.getLogger("polybuilder")


def _configure_logging(verbosity: int) -> None:
    level = logging.WARNING
    if verbosity == 1:
        level = logging.INFO
    elif verbosity >= 2:
        level = logging.DEBUG
    logging.basicConfig(
        level=level,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )


def _add_common_output_args(sp: argparse.ArgumentParser) -> None:
    sp.add_argument(
        "-o", "--output-dir",
        default=".",
        help="Directory for output files (default: current directory).",
    )
    sp.add_argument(
        "--no-editconf",
        action="store_true",
        help="Skip the 'gmx editconf' post-processing step.",
    )
    sp.add_argument(
        "--keep-intermediates",
        action="store_true",
        help="Keep the raw polymer.pdb / polymer.rtp / rearranged_polymer.pdb intermediates.",
    )


def _library_arg(sp: argparse.ArgumentParser) -> None:
    sp.add_argument(
        "--library",
        action="append",
        default=[],
        metavar="FILE",
        help=(
            "Path to a user Python file that registers additional monomers. "
            "May be given multiple times. See README for the format."
        ),
    )


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="polybuilder",
        description="Build 3D polymer structures + GROMACS RTP files.",
    )
    parser.add_argument(
        "-v", "--verbose",
        action="count",
        default=0,
        help="Increase log verbosity (-v INFO, -vv DEBUG).",
    )

    sub = parser.add_subparsers(dest="command", required=True)

    # ---------------- build ----------------
    build_parser = sub.add_parser("build", help="Build a polymer non-interactively.")
    _library_arg(build_parser)
    build_parser.add_argument(
        "--first-residue",
        dest="first_residue",
        default="MMA",
        help=(
            "First (majority) residue in each repeat block. Defaults to MMA. "
            "Any registered residue works (e.g. NIPAM, DABCO)."
        ),
    )
    build_parser.add_argument(
        "--comonomer", "--betaine",
        dest="comonomer",
        required=True,
        help=(
            "Second residue in each repeat block. Built-ins: DMAPS, A3316, MMA. "
            "Additional residues available via --library. "
            "'--betaine' is a legacy alias."
        ),
    )
    build_parser.add_argument(
        "--n-first", "--n-mma",
        dest="n_first", type=int, required=True,
        help=(
            "Number of first-residue units per repeat block "
            "('--n-mma' is a legacy alias, meaningful only when --first-residue=MMA)."
        ),
    )
    build_parser.add_argument(
        "--n-comonomer", "--n-betaine",
        dest="n_comonomer", type=int, required=True,
        help="Number of comonomer units per repeat block ('--n-betaine' is a legacy alias).",
    )
    build_parser.add_argument(
        "--repeats", type=int, required=True,
        help="How many times the (first_n + comonomer_m) block repeats.",
    )
    build_parser.add_argument(
        "--bridge", type=int, default=0,
        help=(
            "Extra bridge methylenes in a sulfobetaine (default: 0). "
            "Must be 0 unless at least one residue is a sulfobetaine."
        ),
    )
    build_parser.add_argument(
        "--cap", action="store_true",
        help="Append one extra first-residue unit at the tail as a cap.",
    )
    build_parser.add_argument(
        "--initiator",
        default=None,
        help=(
            "Optional radical initiator name. Built-ins: KPS, AIBA, TBHP. "
            "Additional initiators available via --library."
        ),
    )
    build_parser.add_argument(
        "--ends",
        default="both",
        choices=("both", "head", "tail", "none"),
        help=(
            "Which chain end(s) to cap with the initiator (default: both). "
            "Ignored if --initiator-at or --initiator-percent is given."
        ),
    )
    build_parser.add_argument(
        "--initiator-at",
        dest="initiator_at",
        default=None,
        metavar="LIST",
        help=(
            "Explicit initiator placement. Comma-separated list of positions: "
            "'head', 'tail', or a 0-based residue index (e.g. 'head,tail,5,12'). "
            "Overrides --ends and --initiator-percent."
        ),
    )
    build_parser.add_argument(
        "--initiator-percent",
        dest="initiator_percent",
        type=float,
        default=None,
        metavar="PCT",
        help=(
            "Target initiator loading as a percentage of the number of "
            "residues. E.g. 5 → ~5%% of residues get an initiator (placed at "
            "head, tail, then evenly-spaced mid-chain). Overrides --ends."
        ),
    )
    _add_common_output_args(build_parser)

    # ---------------- interactive ----------------
    interactive_parser = sub.add_parser(
        "interactive", help="Prompt for parameters like the legacy run_pipeline.py."
    )
    _library_arg(interactive_parser)
    _add_common_output_args(interactive_parser)

    # ---------------- pubchem ----------------
    pub_parser = sub.add_parser(
        "pubchem",
        help="Print a starter config entry for a PubChem SMILES.",
    )
    pub_parser.add_argument("name", help="Residue name to use.")
    pub_parser.add_argument("smiles", help="Vinyl SMILES from PubChem.")

    # ---------------- list-monomers ----------------
    list_parser = sub.add_parser(
        "list-monomers",
        help="Print all registered monomers (after loading --library).",
    )
    _library_arg(list_parser)

    # ---------------- list-initiators ----------------
    li_parser = sub.add_parser(
        "list-initiators",
        help="Print all registered initiators (after loading --library).",
    )
    _library_arg(li_parser)

    # ---------------- inspect ----------------
    ins_parser = sub.add_parser(
        "inspect",
        help="Read a produced PDB and dump per-residue atoms + connectivity.",
    )
    ins_parser.add_argument("pdb", help="Path to the PDB file to inspect.")
    ins_parser.add_argument(
        "-v", "--full",
        dest="inspect_verbose", action="store_true",
        help="Also list every intra-residue neighbour of every heavy atom.",
    )
    ins_parser.add_argument(
        "--residue",
        dest="inspect_residue", default=None,
        help="Only show residues whose name matches (e.g. VBL, NIF).",
    )

    return parser


def _prompt_int(msg: str, default: int | None = None) -> int:
    while True:
        raw = input(msg).strip()
        if not raw and default is not None:
            return default
        try:
            return int(raw)
        except ValueError:
            print("Please enter an integer.", file=sys.stderr)


def _apply_library(args: argparse.Namespace) -> None:
    for path in getattr(args, "library", None) or []:
        load_user_library(path)


def _pick_residue(prompt: str, exclude: str | None = None) -> str:
    residues = [r for r in available_residues() if r != exclude]
    print(f"{prompt}")
    for i, name in enumerate(residues, 1):
        print(f"  {i}. {name}")
    while True:
        raw = input(f"Enter choice (1-{len(residues)}): ").strip()
        try:
            choice = int(raw)
            if 1 <= choice <= len(residues):
                return residues[choice - 1]
        except ValueError:
            pass
        print("Invalid choice; try again.", file=sys.stderr)


def _run_interactive(args: argparse.Namespace) -> int:
    _apply_library(args)
    first = _pick_residue("Select first (majority) residue:")
    comonomer = _pick_residue("Select second residue:", exclude=first)

    n_first = _prompt_int(f"{first} units per block: ")
    n_c = _prompt_int(f"{comonomer} units per block: ")
    reps = _prompt_int("Block repeats: ")
    bridge = _prompt_int("Extra bridge methylenes (0 unless a sulfobetaine is used): ")
    cap = input(f"Add {first} cap? (y/n): ").strip().lower() == "y"

    spec = PolymerSpec(
        first_residue=first,
        comonomer=comonomer,
        n_first=n_first,
        n_comonomer=n_c,
        repeats=reps,
        extra_bridge=bridge,
        cap=cap,
    )

    result = run_pipeline(
        spec,
        output_dir=args.output_dir,
        run_editconf=not args.no_editconf,
        keep_intermediates=args.keep_intermediates,
    )
    print(f"Final PDB: {result.pdb_path}")
    print(f"Final RTP: {result.rtp_path}")
    return 0


def _parse_initiator_at(value: str | None) -> list[str | int] | None:
    """Parse '--initiator-at head,tail,5,12' → ['head','tail',5,12]."""
    if value is None:
        return None
    positions: list[str | int] = []
    for token in value.split(","):
        tok = token.strip()
        if not tok:
            continue
        if tok in ("head", "tail"):
            positions.append(tok)
        else:
            try:
                positions.append(int(tok))
            except ValueError as exc:
                raise PolybuilderError(
                    f"--initiator-at entry {tok!r} is not 'head', 'tail', "
                    f"or an integer"
                ) from exc
    return positions


def _run_build(args: argparse.Namespace) -> int:
    _apply_library(args)
    spec = PolymerSpec(
        first_residue=args.first_residue,
        comonomer=args.comonomer,
        n_first=args.n_first,
        n_comonomer=args.n_comonomer,
        repeats=args.repeats,
        extra_bridge=args.bridge,
        cap=args.cap,
        initiator=args.initiator,
        ends=args.ends,
        initiator_positions=_parse_initiator_at(args.initiator_at),
        initiator_percent=args.initiator_percent,
    )
    result = run_pipeline(
        spec,
        output_dir=args.output_dir,
        run_editconf=not args.no_editconf,
        keep_intermediates=args.keep_intermediates,
    )
    print(f"Final PDB: {result.pdb_path}")
    print(f"Final RTP: {result.rtp_path}")
    return 0


def _run_pubchem(args: argparse.Namespace) -> int:
    smiles, atoms = generate_config_entry(args.name, args.smiles)
    if smiles is None:
        print(f"No polymerisable vinyl bond found in {args.smiles!r}.",
              file=sys.stderr)
        return 1
    print(f"# {args.name} — starter config")
    print(f"smiles = {smiles!r}")
    print("atoms = [")
    for atom in atoms or []:
        print(
            f"    ({atom.name!r}, {atom.ff_type!r}, "
            f"{atom.charge:.4f}, {atom.cgnr}),"
        )
    print("]")
    return 0


def _run_list_monomers(args: argparse.Namespace) -> int:
    _apply_library(args)
    for name in available_residues():
        print(name)
    return 0


def _run_list_initiators(args: argparse.Namespace) -> int:
    _apply_library(args)
    for name in available_initiators():
        print(name)
    return 0


def _run_inspect(args: argparse.Namespace) -> int:
    inspect_pdb_cmd(
        args.pdb,
        verbose=args.inspect_verbose,
        residue_filter=args.inspect_residue,
    )
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    _configure_logging(args.verbose)

    try:
        if args.command == "build":
            return _run_build(args)
        if args.command == "interactive":
            return _run_interactive(args)
        if args.command == "pubchem":
            return _run_pubchem(args)
        if args.command == "list-monomers":
            return _run_list_monomers(args)
        if args.command == "list-initiators":
            return _run_list_initiators(args)
        if args.command == "inspect":
            return _run_inspect(args)
    except PolybuilderError as exc:
        log.error("polybuilder: %s", exc)
        return 2
    except KeyboardInterrupt:
        log.error("Interrupted by user.")
        return 130

    parser.error(f"Unknown command: {args.command}")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
