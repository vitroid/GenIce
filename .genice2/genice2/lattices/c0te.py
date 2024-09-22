#!/usr/bin/python
"""
Usage: genice2 c0te
"""

import genice2.lattices
from genice2.cell import cellvectors
from logging import getLogger
import numpy as np


def usage():
    logger = getLogger()
    logger.info(__doc__)


desc = {
    "ref": {"C0": "Teeratchanan 2015"},
    "usage": usage(),
    "brief": "Filled ice C0 by Teeratchanan (Hydrogen-disordered.) (Positions of guests are supplied.)",
    "test": ({"options": "-r 2 2 2"},),
}


class Lattice(genice2.lattices.Lattice):
    def __init__(self):
        logger = getLogger()

        # Ref. 2atom
        atoms = """
        O1 0.2342 0.4721 0.8019
        O2 0.7648 0.5306 0.2941
        Ne1 -0.0647 0.7868 0.7669
        """

        # Ref.
        # space group: P3_2
        symops = """
          x,            y,            z
         -y,x-y,z+2/3
         -x+y,-x,z+1/3
        """.replace(
            ",", " "
        )

        # Ref. 2cell
        a = 6.177 / 10.0  # nm
        c = 6.054 / 10.0  # nm
        C = 120.0

        self.cell = cellvectors(a, a, c, C=C)

        # helper routines to make from CIF-like data
        from genice2 import CIF

        atomd = CIF.atomdic(atoms)
        atoms = CIF.fullatoms(atomd, CIF.symmetry_operators(symops))

        self.cagetype = []
        self.cagepos = []
        for atomname, pos in atoms:
            if atomname == "Ne1":
                self.cagetype.append(atomname)
                self.cagepos.append(pos)

        self.waters, self.pairs = CIF.waters_and_pairs(
            self.cell, atomd, CIF.symmetry_operators(symops)
        )

        self.density = (
            18 * len(self.waters) / 6.022e23 / (np.linalg.det(self.cell) * 1e-21)
        )
        self.coord = "relative"
