#!/usr/bin/python

import numpy as np
from genice2 import CIF
from genice2.cell import cellvectors
import genice2.lattices
desc = {
    "ref": {
        "engel17": "Engel 2018",
        "DDR": "IZA Database"
    },
    "usage": "No options available.",
    "brief": "Hypothetical zeolitic ice"
}


class Lattice(genice2.lattices.Lattice):
    def __init__(self):
        atoms = """
    O1       0.7281    0.0539    0.0698
    O2       0.1309    0.2618    0.1075
    O3       0.1999    0.3999    0.1718
    O4       0.1187    0.2374    0.2309
    O5       0.2256    0.0000    0.0000
    O6       0.0000    0.0000    0.1956
    O7       0.0000    0.0000    0.1166
        """

        symops = """
+x,+y,+z
2/3+x,1/3+y,1/3+z
1/3+x,2/3+y,2/3+z
-y,+x-y,+z
2/3-y,1/3+x-y,1/3+z
1/3-y,2/3+x-y,2/3+z
-x+y,-x,+z
2/3-x+y,1/3-x,1/3+z
1/3-x+y,2/3-x,2/3+z
-y,-x,+z
2/3-y,1/3-x,1/3+z
1/3-y,2/3-x,2/3+z
-x+y,+y,+z
2/3-x+y,1/3+y,1/3+z
1/3-x+y,2/3+y,2/3+z
+x,+x-y,+z
2/3+x,1/3+x-y,1/3+z
1/3+x,2/3+x-y,2/3+z
-x,-y,-z
2/3-x,1/3-y,1/3-z
1/3-x,2/3-y,2/3-z
+y,-x+y,-z
2/3+y,1/3-x+y,1/3-z
1/3+y,2/3-x+y,2/3-z
+x-y,+x,-z
2/3+x-y,1/3+x,1/3-z
1/3+x-y,2/3+x,2/3-z
+y,+x,-z
2/3+y,1/3+x,1/3-z
1/3+y,2/3+x,2/3-z
+x-y,-y,-z
2/3+x-y,1/3-y,1/3-z
1/3+x-y,2/3-y,2/3-z
-x,-x+y,-z
2/3-x,1/3-x+y,1/3-z
1/3-x,2/3-x+y,2/3-z
        """.replace(',', ' ')

        a = 13.7950 / 10.0  # nm
        b = 13.7950 / 10.0  # nm
        c = 40.7500 / 10.0  # nm
        A = 90
        B = 90
        C = 120

        self.cell = cellvectors(a, b, c, A, B, C)

        # helper routines to make from CIF-like data
        atomd = CIF.atomdic(atoms)
        sops = CIF.symmetry_operators(symops)
        self.waters, self.fixed = CIF.waters_and_pairs(self.cell, atomd, sops)

        # set self.pairs in this way for hydrogen-ordered ices.
        self.pairs = self.fixed

        self.density = 18 * len(self.waters) / 6.022e23 / \
            (np.linalg.det(self.cell) * 1e-21)

        self.coord = "relative"
