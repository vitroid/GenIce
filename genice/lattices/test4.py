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

bondlen = 3

celltype = 'rect'

cell = """
12.747893943706936 12.747893943706936 12.747893943706936
"""

from genice import libgenice, FrankKasper
import logging
logger = logging.getLogger()

cagepos, cagetype = libgenice.parse_cages(cages)
c                 = libgenice.parse_cell(cell, celltype)
waters = [w for w in FrankKasper.toWater(cagepos, c)]
coord = "relative"
density = FrankKasper.estimate_density(waters, c, 2.76)
