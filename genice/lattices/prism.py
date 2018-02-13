# coding: utf-8
"""
Note: Due to the technical limitation in the GenIce algorithm, the minimum lattice size is larger than the crystallographic unit cell size.
"""

from math import sin, pi, cos
import logging

def usage():
    logger = logging.getLogger()
    logger.info("** prism module **")
    logger.info("prism module accepts two arguments.")
    logger.info("prism[sides,rows]")
    logger.info("By default, sides are 6 and rows are 10.")
    logger.info("prism[5] prepares pentagonal prism of 10 rows.")
    logger.info("prism[8,5] prepares octagonal prism of 5 rows.")
    logger.info("------------------")
    


def argparser(arg):
    global sides, rows, bondlen, density, cell, celltype, waters, coord
    a = [int(x) for x in arg.split(",")]
    if len(a)>0:
        sides = a[0]
    if len(a)>1:
        rows = a[1]
    L = 2.75
    bondlen = 3
    R = L/2/sin(pi/sides)
    celltype = "rect"
    cell = [L*20, L*20, L*rows]
    density = sides*rows / (L**3 * 400 * rows) * 18 / 6.022e23 * 1e24

    waters = []
    for j in range(rows):
        for i in range(sides):
            x = R * cos(i*pi*2/sides)
            y = R * sin(i*pi*2/sides)
            z = j * L
            waters.append([x,y,z])

    coord = "absolute"


# default.
argparser("6,10")
