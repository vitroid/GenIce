# coding: utf-8
"""
Test for generating water positions from cage positions
"""


cages="""
12 0.5 0.5 0.5
14 0.5 0.0 -0.25
14 0.0 0.25 0.5
14 0.75 0.5 0.0
14 0.5 0.0 0.25
12 0.0 0.0 0.0
14 0.25 0.5 0.0
14 0.0 -0.25 0.5
"""



from genice         import FrankKasper
from genice.lattice import parse_cages
from genice.cell import cellvectors
cell = cellvectors(a=12.747893943706936,
                   b=12.747893943706936,
                   c=12.747893943706936)
cagepos, cagetype = parse_cages(cages)
waters = [w for w in FrankKasper.toWater(cagepos, cell)]
coord = "relative"
density = FrankKasper.estimate_density(waters, cell, 2.76)
bondlen = 2.76 * 1.2

