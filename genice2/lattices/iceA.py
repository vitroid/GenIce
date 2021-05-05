#!/usr/bin/python
# coding: utf-8

import numpy as np
from genice2 import CIF
from genice2.cell import cellvectors
import genice2.lattices
desc = {"ref": {"A": 'Baez 1998'},
        "usage": "No options available.",
        "brief": "Hypothetical ice A."
        }


class Lattice(genice2.lattices.Lattice):
    def __init__(self):
        atoms = """
        O1 0.0995 0.1901 0.2835
        O2 0.3431 0.3469 0.7598
        O3 0.4055 0.0558 0.5013
        H1 -0.0149 0.1524 0.3658
        H2 0.0921 0.3365 0.2570
        H3 0.4441 0.3484 0.8630
        H4 0.3713 0.2326 0.6746
        H5 0.3045 0.1078 0.4119
        H6 0.3847 -0.0897 0.5194
        """

        # P4_1, No. 76
        symops = """
              x            y            z
             -x           -y          1/2+z
             -y            x          1/4+z
              y           -x          3/4+z
        """

        # in nm
        a, b, c = 0.6733, 0.6733, 0.7164

        self.cell = cellvectors(a, b, c)

        # helper routines to make from CIF-like data
        atomd = CIF.atomdic(atoms)
        sops = CIF.symmetry_operators(symops)
        self.waters, self.fixed = CIF.waters_and_pairs(self.cell, atomd, sops)

        # set self.pairs in this way for hydrogen-ordered ices.
        self.pairs = self.fixed

        self.density = 18 * len(self.waters) / 6.022e23 / \
            (np.linalg.det(self.cell) * 1e-21)

        self.coord = "relative"
