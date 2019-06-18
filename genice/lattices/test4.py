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
bondlen = 2.76*1.1

