# coding: utf-8
"""
Test for generating water positions from cage positions
"""


import genice.lattices
from genice         import FrankKasper
from genice.genice  import parse_cages
from genice.cell import cellvectors

class Lattice(genice.lattices.Lattice):
    def __init__(self):
        self.cages = """
        12    0.5000    0.2500    0.2500
        12    0.5000    0.5000    0.5000
        12    0.2500    0.5000    0.2500
        12    0.0000    0.2500    0.7500
        12    0.2500    0.0000    0.7500
        12    0.0000    0.5000    0.0000
        12    0.2500    0.2500    0.5000
        12    0.5000    0.0000    0.0000
        12    0.7500    0.5000    0.7500
        12    0.0000    0.0000    0.5000
        12    0.0000    0.7500    0.2500
        12    0.2500    0.7500    0.0000
        12    0.7500    0.2500    0.0000
        12    0.7500    0.0000    0.2500
        12    0.5000    0.7500    0.7500
        12    0.7500    0.7500    0.5000
        16    0.6250    0.1250    0.6250
        16    0.1250    0.6250    0.6250
        16    0.8750    0.3750    0.3750
        16    0.3750    0.8750    0.3750
        16    0.3750    0.3750    0.8750
        16    0.6250    0.6250    0.1250
        16    0.8750    0.8750    0.8750
        16    0.1250    0.1250    0.1250
        """



        self.cell = cellvectors(a=12.747893943706936,
                           b=12.747893943706936,
                           c=12.747893943706936)
        cagepos, cagetype = parse_cages(self.cages)
        self.waters = [w for w in FrankKasper.toWater(cagepos, self.cell)]
        self.coord = "relative"
        self.density = FrankKasper.estimate_density(self.waters, self.cell, 2.76)
        self.bondlen = 2.76*1.1

