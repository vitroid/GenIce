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
