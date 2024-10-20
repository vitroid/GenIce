#!/usr/bin/python
# coding: utf-8

import numpy as np
from genice2 import CIF
from genice2.cell import cellvectors
import genice2.lattices
from logging import getLogger

desc = {
    "ref": {"M": "Mochizuki 2024"},
    "usage": "No options available.",
    "brief": "A hypothetical hydrogen-ordered high-density ice.",
    "test": ({"options": "--depol=none"},),
}


class Lattice(genice2.lattices.Lattice):
    def __init__(self):
        logger = getLogger()

        # only the hydrogen atoms whose locations are confirmed are listed.
        atoms = """
O1     0.34919  0.81250  0.29286  1.00000 Dx,Dy,Dz 
H1     0.51098  0.72368  0.32678  1.00000 Dx,Dy,Dz 
H2     0.45942 -0.08553  0.25446  1.00000 Dx,Dy,Dz 
O2     0.83457  0.87829  0.89464  1.00000 Dx,Dy,Dz 
H3     0.75474 -0.02895 -0.02053  1.00000 Dx,Dy,Dz 
H4     0.00714  0.82829  0.00893  1.00000 Dx,Dy,Dz 
O3     0.37618  0.61447  0.80625  1.00000 Dx,Dy,Dz 
H5     0.21915  0.63619  0.66429  1.00000 Dx,Dy,Dz 
H6     0.52378  0.70921  0.81429  1.00000 Dx,Dy,Dz 
O4     0.14808  0.05789  0.58482  1.00000 Dx,Dy,Dz 
H7     0.03727  0.00066  0.69196  1.00000 Dx,Dy,Dz 
H8     0.21503 -0.03684  0.49553  1.00000 Dx,Dy,Dz 
        """

        # space group: P 1 21 1 No. 4
        # http://img.chem.ucl.ac.uk/sgp/large/004ay1.htm
        symops = """
x,y,z
-x,y+1/2,-z
    """.translate(
            {ord(","): " "}
        )

        a = 4.3 / 10.0  # nm
        b = 7.6 / 10.0  # nm
        c = 5.727 / 10.0  # nm
        A = 90
        B = 102.09
        C = 90

        self.cell = cellvectors(a, b, c, A, B, C)

        # helper routines to make from CIF-like data
        atomd = CIF.atomdic(atoms)
        sops = CIF.symmetry_operators(symops)
        self.waters, self.fixed, self.pairs = CIF.waters_and_pairs(
            self.cell, atomd, sops, rep=[2, 1, 2], partial_order=True
        )

        self.density = (
            18 * len(self.waters) / 6.022e23 / (np.linalg.det(self.cell) * 1e-21)
        )
        self.bondlen = 0.3

        self.coord = "relative"
