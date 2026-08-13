"""Typed exceptions raised by :mod:`polybuilder`."""

from __future__ import annotations


class PolybuilderError(Exception):
    """Base class for all polybuilder errors."""


class InvalidSmilesError(PolybuilderError):
    """Raised when RDKit cannot parse a SMILES string."""


class UnknownMonomerError(PolybuilderError):
    """Raised when a requested monomer name is not in the library."""


class InvalidPolymerSpecError(PolybuilderError):
    """Raised when a :class:`PolymerSpec` has non-sensical values."""


class EmbeddingError(PolybuilderError):
    """Raised when RDKit fails to embed a 3D conformer."""


class GromacsNotFoundError(PolybuilderError):
    """Raised when the ``gmx`` executable is required but not on PATH."""
