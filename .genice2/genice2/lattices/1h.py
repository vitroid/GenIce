# coding: utf-8

from genice2.cell import cellvectors
import genice2.lattices
desc = {"ref": {},
        "usage": "No options available.",
        "brief": "Most popular Ice I (hexagonal)"
        }


class Lattice(genice2.lattices.Lattice):
    def __init__(self):
        self.density = 0.92  # default self.density

        self.bondlen = 3  # bond threshold

        self.cell = cellvectors(
            a=7.84813412606925, b=7.37735062301457, c=9.06573834219084)

        self.waters = """
        1.328 1.802 3.38
        5.267 4.524 1.109
        6.58 5.442 3.365
        5.267 4.542 5.629
        2.623 0.877 5.644
        2.667 5.488 5.625
        5.241 1.756 1.12
        5.241 1.774 5.64
        1.354 4.588 7.888
        1.354 4.57 3.369
        2.667 5.47 1.105
        2.623 0.858 1.124
        6.537 0.831 3.384
        6.537 0.849 7.903
        6.581 5.461 7.884
        1.328 1.82 7.899
        """

        self.coord = "absolute"
