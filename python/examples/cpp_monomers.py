"""Example plugin library — now shipped as built-ins.

Everything that used to be defined here (A3367, A3361, M3295, NIPAM, IBOA,
VBD/DABCO, BVD and the APS initiator) is now a **built-in** monomer, so you no
longer need ``--library`` to use them::

    polybuilder homo DABCO -n 20
    polybuilder copo DMAPS A3367 -n 10

This file is kept as a working demonstration of the plugin mechanism and for
backward compatibility with older scripts that do::

    polybuilder build --library examples/cpp_monomers.py --comonomer A3367 ...

It simply re-exports the built-in definitions.  Loading it is now a harmless
no-op.  To add a genuinely new monomer, define your own ``ResidueFragment`` /
``BetaineBase`` and put it in the ``MONOMERS`` dict below — see the README.
"""

from polybuilder.library.catalog import (
    A3361_BASE,
    A3367_BASE,
    APS_INITIATOR,
    BVD_RESIDUE,
    IBOA_RESIDUE,
    M3295_BASE,
    NIPAM_RESIDUE,
    VBD_RESIDUE,
)

MONOMERS = {
    # sulfobetaines
    "A3367": A3367_BASE,
    "M3295": M3295_BASE,
    "A3361": A3361_BASE,
    # single-fragment monomers (extra_bridge must be 0)
    "NIPAM": NIPAM_RESIDUE,
    "IBOA": IBOA_RESIDUE,
    "VBD": VBD_RESIDUE,
    # legacy alias for the vinyl-benzyl-DABCO monomer
    "DABCO": VBD_RESIDUE,
    "BVD": BVD_RESIDUE,
}

INITIATORS = {
    "APS": APS_INITIATOR,
}
