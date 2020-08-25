# coding: utf-8
"""
Test for generating water positions from cage positions
"""


import genice.lattices
from genice         import FrankKasper
from genice.valueparsers import parse_cages
from genice.cell    import cellvectors

class Lattice(genice.lattices.Lattice):
    def __init__(self):
        self.cages="""
        12 0.5 0.5 0.5
        14 0.5 0.0 -0.25
        14 0.0 0.25 0.5
        14 0.75 0.5 0.0
        14 0.5 0.0 0.25
        12 0.0 0.0 0.0
        14 0.25 0.5 0.0
        14 0.0 -0.25 0.5
        """



        self.cell = cellvectors(a=12.747893943706936,
                           b=12.747893943706936,
                           c=12.747893943706936)
        cagepos, cagetype = parse_cages(self.cages)
        self.waters = [w for w in FrankKasper.toWater(cagepos, self.cell)]
        self.coord = "relative"
        self.density = FrankKasper.estimate_density(self.waters, self.cell, 2.76)
        self.bondlen = 2.76 * 1.2

