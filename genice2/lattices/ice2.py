#!/usr/bin/python
"""
Usage: genice2 two
"""
from logging import getLogger
from genice2 import CIF
import numpy as np
from genice2.cell import cellvectors
import genice2.lattices


def usage():
    logger = getLogger()
    logger.info(__doc__)


desc = {"ref": {"2atom": "Kamb 2003",
                "2cell": "Kamb 1964",
                "C1": "Londono 1988"},
        "usage": usage(),
        "brief": "Hydrogen-ordered ice II."
        }


class Lattice(genice2.lattices.Lattice):
    def __init__(self):
        logger = getLogger()

        # Ref. 2atom
        atoms = """
        O1   0.2716    0.0259   -0.1471
        O2   0.4798    0.7571    0.3389
        D1   0.7284    0.4038    0.4034
        D2   0.1491    0.0412   -0.2023
        D3   0.7420    0.1978    0.3708
        D4   0.4232    0.1954   -0.0164
        """

        # Ref. 2cell
        # space group: R-3
        symops = """
          x,            y,            z
          z,            x,            y
          y,            z,            x

         -x,           -y,           -z
         -z,           -x,           -y
         -y,           -z,           -x
        """.replace(',', ' ')

        # Ref. 2cell
        a = 7.78 / 10.0  # nm
        A = 113.1

        self.cell = cellvectors(a, a, a, A, A, A)

        # helper routines to make from CIF-like data
        atomd = CIF.atomdic(atoms)
        sops = CIF.symmetry_operators(symops)
        self.waters, self.fixed = CIF.waters_and_pairs(self.cell, atomd, sops)

        # set pairs in this way for hydrogen-ordered ices.
        self.pairs = self.fixed

        self.density = 18 * len(self.waters) / 6.022e23 / \
            (np.linalg.det(self.cell) * 1e-21)

        self.coord = "relative"
