#!/usr/bin/python
# coding: utf-8

import numpy as np
from genice2 import CIF
from genice2.cell import cellvectors
import genice2.lattices
from logging import getLogger

desc = {
    "ref": {"hs6": "Matsumoto 2021"},
    "usage": "No options available.",
    "brief": "A ultralow-density ice with channels. It is actually an AFI.",
    "test": ({"options": "--depol=none"},),
}


class Lattice(genice2.lattices.Lattice):
    def __init__(self):
        logger = getLogger()

        # The locations are estimated by hand.
        atoms = """
O1     0.44444  0.11111  0.2      
        """

        # space group: P 6/m cc No. 192
        # http://img.chem.ucl.ac.uk/sgp/large/192az3.htm
        symops = """
     x,            y,            z
     -y,           x-y,           z
    -x+y,          -x,            z
     -x,           -y,            z
     x-y,           x,            z
      y,          -x+y,           z
     -y,           -x,          1/2+z
    -x+y,           y,          1/2+z
      x,           x-y,         1/2+z
      y,            x,          1/2+z
     x-y,          -y,          1/2+z
     -x,          -x+y,         1/2+z
 
     -x,           -y,           -z
      y,          -x+y,          -z
     x-y,           x,           -z
      x,            y,           -z
    -x+y,          -x,           -z
     -y,           x-y,          -z
      y,            x,          1/2-z
     x-y,          -y,          1/2-z
     -x,          -x+y,         1/2-z
     -y,           -x,          1/2-z
    -x+y,           y,          1/2-z
      x,           x-y,         1/2-z
""".translate(
            {ord(","): " "}
        )

        a = 4.5 / 10.0  # nm
        c = 2.5 / 10.0  # nm
        A = 90
        C = 120

        self.cell = cellvectors(a, a, c, A, A, C)

        # helper routines to make from CIF-like data
        atomd = CIF.atomdic(atoms)
        sops = CIF.symmetry_operators(symops)
        self.waters, self.pairs = CIF.waters_and_pairs(self.cell, atomd, sops)

        # self.density = (
        #     18 * len(self.waters) / 6.022e23 / (np.linalg.det(self.cell) * 1e-21)
        # )
        # self.density = 0.9
        self.bondlen = 0.11

        self.coord = "relative"
