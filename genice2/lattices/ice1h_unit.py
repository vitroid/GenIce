# coding: utf-8
#This is not useful for generating the graph

import genice2.lattices
from genice2.cell import cellvectors

class Lattice(genice2.lattices.Lattice):
    def __init__(self):
        self.density = 0.92     #default self.density

        self.bondlen = 3        #bond threshold
        self.cell = """
        4.5328691711 7.84813412606925 7.37735062301457
        """

        self.waters = """
        0.25 0.333 0.125
        0.75 0.833 0.125
        0.25 0.667 0.25
        0.75 0.167 0.25
        0.25 0.667 0.625
        0.75 0.167 0.625
        0.75 0.833 0.75
        0.25 0.333 0.75
        """

        self.coord = "relative"

        self.cell = cellvectors(a=4.5328691711,
                           b=7.84813412606925,
                           c=7.37735062301457)
