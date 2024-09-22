#!/usr/bin/python
# coding: utf-8

import numpy as np
from genice2 import CIF
from genice2.cell import cellvectors
import genice2.lattices
from logging import getLogger

desc = {
    "ref": {"XIV": "Salzmann 2006"},
    "usage": "No options available.",
    "brief": "Ice XIV, a partially hydrogen-ordered counterpart of ice XII. Note that it does not reproduce the occupancies (probability of occupation) of the possible hydrogen sites.",
    "test": ({"options": "--depol=none"},),
}


class Lattice(genice2.lattices.Lattice):
    def __init__(self):
        logger = getLogger()

        # From Table 2 of Salzmann 2006

        # only the hydrogen atoms whose locations are confirmed are listed.
        atoms = """
O1 0.0059(3) 0.2568(5) 0.1304(7) 1.53(2) 1.0000 
O2 0.6308(3) -0.0078(3) 0.2485(7) 1.53(2) 1.0000 
O3 0.2525(4) 0.8858(3) 0.0063(6) 1.53(2) 1.0000 
D6 0.0920(2) 0.2056(3) 0.2671(6) 2.03(1) 1.0000
D9 0.7895(3) 0.9679(3) 0.8954(7) 2.03(1) 1.0000
D11 0.7340(3) 0.4630(3) 0.3225(6) 2.03(1) 1.0000
D15 0.8472(3) 0.3248(3) 0.4010(6) 2.03(1) 1.0000
        """
        # 以下の辺は不確定。
        # D4 0.0557(7) 0.3284(7) 0.9845(1) 2.03(1) 0.407(3)
        # D5 0.5275(5) 0.8410(4) 0.4684(1) 2.03(1) 0.620(4)
        # D12 0.4111(4) 0.5790(5) 0.3625(1) 2.03(1) 0.593(3)
        # D13 0.9018(8) 0.1018(7) 0.8552(2) 2.03(1) 0.380(4)

        # space group: P2_1 2_1 2_1 No. 19 (b)
        # http://img.chem.ucl.ac.uk/sgp/large/019bz1.htm
        symops = """
          x,            y,            z
        1/2+x,        1/2-y,         -z
         -x,          1/2+y,        1/2-z
        1/2-x,         -y,          1/2+z
    """.translate(
            {ord(","): ""}
        )

        a = 8.3499 / 10.0  # nm
        b = 8.1391 / 10.0  # nm
        c = 4.0825 / 10.0  # nm
        A = 90
        B = 90
        C = 90

        self.cell = cellvectors(a, b, c, A, B, C)

        # helper routines to make from CIF-like data
        atomd = CIF.atomdic(atoms)
        sops = CIF.symmetry_operators(symops)
        self.waters, self.fixed, self.pairs = CIF.waters_and_pairs(
            self.cell, atomd, sops, rep=(1, 1, 2), partial_order=True
        )

        self.density = (
            18 * len(self.waters) / 6.022e23 / (np.linalg.det(self.cell) * 1e-21)
        )
        self.bondlen = 0.3

        self.coord = "relative"
