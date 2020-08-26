# coding: utf-8
"""
[sVII] Jeffrey, G A. “Hydrate Inclusion Compounds.” Inclusion Compounds 1 (1984): 135–190.
[CS4] Kosyakov, Viktor I, and T M Polyanskaya. “Using Structural Data for Estimating the Stability of Water Networks in Clathrate and Semiclathrate Hydrates.” Journal of Structural Chemistry 40.2 (1999): 239–245.
[SOD] http://www.iza-structure.org/databases/
"""

import genice2.lattices
from genice2.cell import cellvectors

class Lattice(genice2.lattices.Lattice):
    def __init__(self):
        self.cell = """
        7.8064588643 7.8064588643 7.8064588643
        """
        self.waters = """
        0.0 5.85484414822 3.90322943215
        3.90322943215 1.95161471607 0.0
        1.95161471607 0.0 3.90322943215
        5.85484414822 3.90322943215 0.0
        0.0 3.90322943215 1.95161471607
        3.90322943215 0.0 5.85484414822
        0.0 3.90322943215 5.85484414822
        3.90322943215 0.0 1.95161471607
        3.90322943215 5.85484414822 0.0
        0.0 1.95161471607 3.90322943215
        5.85484414822 0.0 3.90322943215
        1.95161471607 3.90322943215 0.0
        """
        self.coord = "absolute"
        self.bondlen = 3
        self.density = 0.75396428378

        self.cages="""
        K 0 0 0
        K 0.5 0.5 0.5
        """

        self.cell = cellvectors(a=7.8064588643,
                           b=7.8064588643,
                           c=7.8064588643)
