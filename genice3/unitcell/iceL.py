#!/usr/bin/python
"""
Usage: genice2 c0te
"""

import genice3.unitcell
from genice3.util import (
    cellvectors,
    atomdic,
    fullatoms,
    symmetry_operators,
    waters_and_pairs,
    density_in_g_cm3,
)
from logging import getLogger
import numpy as np


desc = {
    "ref": {"L": "Lei 2025"},
    "usage": __doc__,
    "brief": "The hypothetical Ice L",
    "test": ({"options": "-r 2 2 2"},),
}


class UnitCell(genice3.unitcell.UnitCell):
    def __init__(self):

        atoms = """
        O1 0.3488 0.9861 0.3123
        O2 0.9710 0.6650 0.0774
        O3 0.3728 0.7575 0.1387
        """

        # space group: Pbca
        symops = """
        x,            y,            z
        1/2-x,        1/2+y,          z
        x,          1/2-y,        1/2+z
        1/2+x,          y,          1/2-z
    
        -x,           -y,           -z
        1/2+x,        1/2-y,         -z
        -x,          1/2+y,        1/2-z
        1/2-x,         -y,          1/2+z
        """.replace(
            ",", " "
        )

        # Ref. 2cell
        a = 6.23 / 10  # nm
        b = 7.24 / 10  # nm
        c = 11.88 / 10  # nm

        cell = cellvectors(a, b, c)

        atomd = atomdic(atoms)
        waters, _ = waters_and_pairs(cell, atomd, symmetry_operators(symops))

        density = density_in_g_cm3(len(waters), cell)

        coord = "relative"
        super().__init__(
            cell=cell,
            waters=waters,
            density=density,
            coord=coord,
            bondlen=0.3,
        )
