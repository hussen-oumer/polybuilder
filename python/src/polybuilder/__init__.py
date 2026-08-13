"""polybuilder — polymer PDB/RTP builder for GROMACS.

Public API
----------
    PolymerSpec            — parameters for one polymer
    PolymerResult          — result of :func:`build_polymer`
    homopolymer(mono, n)   — quick PolymerSpec for a homopolymer
    copolymer(a, b, n)     — quick PolymerSpec for a two-monomer copolymer
    build_polymer(spec)    — builds one polymer, writes .pdb + .rtp
    run_pipeline(spec)     — full pipeline: build + interleave + reorder + gmx
    PipelineResult         — result of :func:`run_pipeline`
"""

from __future__ import annotations

from .core.builder import (
    PolymerResult,
    PolymerSpec,
    build_polymer,
    copolymer,
    homopolymer,
)
from .core.pipeline import PipelineResult, run_pipeline
from .library import (
    available_betaines,
    available_initiators,
    available_residues,
    get_betaine,
    get_initiator,
    get_residue,
    load_user_library,
    register_betaine,
    register_initiator,
    register_residue,
)
from .library.config import AtomSpec, BetaineBase, Initiator, ResidueFragment
from .support.exceptions import (
    EmbeddingError,
    GromacsNotFoundError,
    InvalidPolymerSpecError,
    InvalidSmilesError,
    PolybuilderError,
    UnknownMonomerError,
)
from .support.inspect import inspect, inspect_pdb, print_report

__version__ = "0.4.0"

__all__ = [
    "PolymerSpec",
    "PolymerResult",
    "build_polymer",
    "homopolymer",
    "copolymer",
    "PipelineResult",
    "run_pipeline",
    "AtomSpec",
    "BetaineBase",
    "Initiator",
    "ResidueFragment",
    "register_residue",
    "get_residue",
    "available_residues",
    "register_initiator",
    "get_initiator",
    "available_initiators",
    "load_user_library",
    "inspect", "inspect_pdb", "print_report",
    # deprecated aliases
    "register_betaine",
    "get_betaine",
    "available_betaines",
    "PolybuilderError",
    "InvalidSmilesError",
    "InvalidPolymerSpecError",
    "UnknownMonomerError",
    "EmbeddingError",
    "GromacsNotFoundError",
    "__version__",
]
