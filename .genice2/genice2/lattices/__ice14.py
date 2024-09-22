#!/usr/bin/python

import numpy as np
from genice2 import CIF
from genice2.cell import cellvectors
import genice2.lattices
desc = {
    "ref": {
        "14": 'Salzmann, C. G., Radaelli, P., Hallbrucker, A. & Mayer, E. The Preparation and Structures of Hydrogen Ordered Phases of Ice. Science 311, 1758–1761 (2006).'},
    "usage": "No options available.",
    "brief": "Ice XIV. (EXPERIMENTAL)"}

assert False, "IT DOES NOT WORK."


class Lattice(genice2.lattices.Lattice):
    def __init__(self):
        atoms = """
        O1 0.0059(3) 0.2568(5) 0.1304(7) 1.53(2) 1.0000
        O2 0.6308(3) -0.0078(3) 0.2485(7) 1.53(2) 1.0000
        O3 0.2525(4) 0.8858(3) 0.0063(6) 1.53(2) 1.0000
        D5 0.5275(5) 0.8410(4) 0.4684(1) 2.03(1) 0.620(4)
        D6 0.0920(2) 0.2056(3) 0.2671(6) 2.03(1) 1.0000
        D9 0.7895(3) 0.9679(3) 0.8954(7) 2.03(1) 1.0000
        D11 0.7340(3) 0.4630(3) 0.3225(6) 2.03(1) 1.0000
        D12 0.4111(4) 0.5790(5) 0.3625(1) 2.03(1) 0.593(3)
        D15 0.8472(3) 0.3248(3) 0.4010(6) 2.03(1) 1.0000
        """
        # Minor hydrogen sites omitted.
        # D4 0.0557(7) 0.3284(7) 0.9845(1) 2.03(1) 0.407(3)
        # D13 0.9018(8) 0.1018(7) 0.8552(2) 2.03(1) 0.380(4)

        # space group: P2_1 2_1 2_1 No. 19
        symops = """
             x,            y,            z
            1/2+x,         -y,          1/2-z
            1/2-x,        1/2+y,         -z
             -x,          1/2-y,        1/2+z
        """
        # http://img.chem.ucl.ac.uk/sgp/large/019az3.htm
        #x,            y,            z
        #    1/2+x,        1/2-y,         -z
        #     -x,          1/2+y,        1/2-z
        #    1/2-x,         -y,          1/2+z

        a = 8.3499 / 10.0  # nm
        b = 8.1391 / 10.0  # nm
        c = 4.0825 / 10.0  # nm

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
