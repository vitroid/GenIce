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



celltype = 'rect'

cell = """
12.747893943706936 12.747893943706936 12.747893943706936
"""

from genice         import FrankKasper
from genice.lattice import parse_cages
from genice.cells   import Cell
cagepos, cagetype = parse_cages(cages)
cell9             = Cell(cell, celltype)
waters = [w for w in FrankKasper.toWater(cagepos, cell9.mat)]
coord = "relative"
density = FrankKasper.estimate_density(waters, cell9.mat, 2.76)
bondlen = 2.76 * 1.2
