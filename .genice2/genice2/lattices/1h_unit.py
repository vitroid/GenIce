# coding: utf-8
import genice2.lattices
from genice2.cell import cellvectors


class Lattice(genice2.lattices.Lattice):
    def __init__(self):
        self.density = 0.92  # default self.density

        self.bondlen = 3  # bond threshold

        self.waters = """
        0.125 0.333 0.125
        0.375 0.833 0.125
        0.125 0.667 0.25
        0.375 0.167 0.25
        0.125 0.667 0.625
        0.375 0.167 0.625
        0.375 0.833 0.75
        0.125 0.333 0.75

        0.625 0.333 0.125
        0.875 0.833 0.125
        0.625 0.667 0.25
        0.875 0.167 0.25
        0.625 0.667 0.625
        0.875 0.167 0.625
        0.875 0.833 0.75
        0.625 0.333 0.75
        """

        self.coord = "relative"

        self.cell = cellvectors(a=4.5328691711 * 2,
                                b=7.84813412606925,
                                c=7.37735062301457)
