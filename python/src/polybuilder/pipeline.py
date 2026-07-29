"""End-to-end orchestration: build → interleave PDB → reorder RTP → gmx editconf."""

from __future__ import annotations

import logging
import os
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

from .builder import PolymerSpec, build_polymer
from .exceptions import GromacsNotFoundError
from .pdb import interleave_pdb
from .rtp import reorder_rtp_by_cgnr

log = logging.getLogger(__name__)


@dataclass
class PipelineResult:
    """Files produced by :func:`run_pipeline`."""

    pdb_path: Path
    rtp_path: Path
    intermediate_files: list


def run_pipeline(
    spec: PolymerSpec,
    output_dir: os.PathLike | str = ".",
    run_editconf: bool = True,
    keep_intermediates: bool = False,
    final_pdb_name: str = "poly.pdb",
    final_rtp_name: str = "rearranged_polymer.rtp",
) -> PipelineResult:
    """Run the full build pipeline into ``output_dir`` and return the paths.

    * Builds a polymer from ``spec``.
    * Interleaves heavy/H atoms in the PDB.
    * Reorders the RTP by charge-group number.
    * Optionally runs ``gmx editconf`` to centre the molecule and put it in a
      periodic box.  Requires ``gmx`` on the PATH.
    """
    out = Path(output_dir).resolve()
    out.mkdir(parents=True, exist_ok=True)

    raw_pdb = out / "polymer.pdb"
    raw_rtp = out / "polymer.rtp"
    inter_pdb = out / "rearranged_polymer.pdb"
    final_rtp = out / final_rtp_name
    final_pdb = out / final_pdb_name

    log.info("Building polymer with %s", spec)
    build_polymer(spec, str(raw_pdb), str(raw_rtp))

    log.info("Interleaving heavy atoms and hydrogens in PDB")
    interleave_pdb(str(raw_pdb), str(inter_pdb))

    log.info("Reordering RTP by charge-group number")
    reorder_rtp_by_cgnr(str(raw_rtp), str(final_rtp))

    intermediates = [raw_pdb, raw_rtp, inter_pdb]

    if run_editconf:
        _run_editconf(inter_pdb, final_pdb)
    else:
        # No gmx step: promote the interleaved PDB to the "final" name.
        shutil.copyfile(inter_pdb, final_pdb)

    if not keep_intermediates:
        for path in intermediates:
            try:
                if path.exists() and path != final_pdb and path != final_rtp:
                    path.unlink()
            except OSError as exc:
                log.warning("Could not remove intermediate %s: %s", path, exc)
        intermediates = []

    return PipelineResult(
        pdb_path=final_pdb,
        rtp_path=final_rtp,
        intermediate_files=intermediates,
    )


def _run_editconf(input_pdb: Path, output_pdb: Path) -> None:
    """Run ``gmx editconf`` to centre and box the molecule."""
    gmx = shutil.which("gmx")
    if gmx is None:
        raise GromacsNotFoundError(
            "gmx executable not found on PATH. Install GROMACS or re-run with "
            "run_editconf=False."
        )

    cmd = [gmx, "editconf", "-f", str(input_pdb), "-o", str(output_pdb)]
    log.info("Running: %s", " ".join(cmd))
    try:
        result = subprocess.run(
            cmd,
            check=True,
            capture_output=True,
            text=True,
        )
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(
            f"gmx editconf failed (exit {exc.returncode}):\n{exc.stderr}"
        ) from exc
    log.debug("gmx editconf stdout: %s", result.stdout)


__all__ = ["PipelineResult", "run_pipeline"]
