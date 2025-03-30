#!/usr/bin/python
"""
Usage: genice2 c0te
"""

import genice2.lattices
from genice2.cell import cellvectors
from logging import getLogger
import numpy as np


desc = {
    "ref": {"L": "Lei 2025"},
    "usage": __doc__,
    "brief": "The hypothetical Ice L",
    "test": ({"options": "-r 2 2 2"},),
}


class Lattice(genice2.lattices.Lattice):
    def __init__(self):
        logger = getLogger()

        # Ref. 2atom
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

        self.cell = cellvectors(a, b, c)

        # helper routines to make from CIF-like data
        from genice2 import CIF

        atomd = CIF.atomdic(atoms)
        atoms = CIF.fullatoms(atomd, CIF.symmetry_operators(symops))

        self.waters, self.pairs = CIF.waters_and_pairs(
            self.cell, atomd, CIF.symmetry_operators(symops)
        )

        self.density = (
            18 * len(self.waters) / 6.022e23 / (np.linalg.det(self.cell) * 1e-21)
        )
        self.coord = "relative"
