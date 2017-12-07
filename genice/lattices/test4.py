# coding: utf-8
"""
Test for generating water positions from cage positions
"""


cages = """
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
bondlen = 2.76*1.1
