#!/usr/bin/python

import numpy as np
from genice2 import CIF
from genice2.cell import cellvectors
import genice2.lattices
desc = {"ref": {"B": 'Baez 1998'},
        "usage": "No options available.",
        "brief": "Hypothetical ice B."
        }


class Lattice(genice2.lattices.Lattice):
    def __init__(self):
        atoms = """
        O1 0.2382 0.7981 0.3422
        O2 0.5000 0.5000 0.0446
        H1 -0.2603 0.9742 0.4613
        H2 0.1446 0.8545 0.2013
        H3 0.4108 0.6151 0.1577
        """

        # P2_1 2_1 2
        symops = """
             x            y            z
            1/2+x        1/2-y         -z
            1/2-x        1/2+y         -z
             -x           -y            z
        """

        # in nm
        a, b, c = 0.7176, 0.4428, 0.5040

        self.cell = cellvectors(a, b, c)

        # helper routines to make from CIF-like data
        atomd = CIF.atomdic(atoms)
        sops = CIF.symmetry_operators(symops)
        # the unit self.cell is too small to handle; multiply (2,2,2)
        self.waters, self.fixed = CIF.waters_and_pairs(
            self.cell, atomd, sops, rep=(2, 2, 2))

        # set self.pairs in this way for hydrogen-ordered ices.
        self.pairs = self.fixed

        self.density = 18 * len(self.waters) / 6.022e23 / \
            (np.linalg.det(self.cell) * 1e-21)

        self.coord = "relative"
