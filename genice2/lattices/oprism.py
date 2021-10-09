# coding: utf-8
"""
Note: Due to the technical limitation in the GenIce algorithm, the minimum lattice size is larger than the crystallographic unit cell size.
"""

import genice2.lattices
from genice2.cell import cellvectors
from logging import getLogger
from math import sin, pi, cos
desc = {
    "ref": {
        "prism": "Koga 2001"
    },
    "usage": "No options available.",
    "brief": "Hydrogen-ordered ice nanotubes."
}


def usage():
    logger = getLogger()
    logger.info("** prism module **")
    logger.info("prism module accepts two arguments.")
    logger.info("prism[sides,rows]")
    logger.info("By default, sides are 6 and rows are 10.")
    logger.info("Number of rows must be even.")
    logger.info("prism[5] prepares pentagonal prism of 10 rows.")
    logger.info("prism[8,6] prepares octagonal prism of 6 rows.")
    logger.info("------------------")


class Lattice(genice2.lattices.Lattice):
    """
Generate a hydrogen-ordered prism ice.

Options:
  rows=r
  sides=s
    """

    def __init__(self, **kwargs):
        logger = getLogger()
        # global sides, rows, bondlen, density, cell, waters, coord
        sides = 6
        rows = 10

        for k, v in kwargs.items():
            if k == "rows":
                rows = int(v)
                assert rows%2 == 0, "Number of rows nums be even."
            elif k == "sides":
                sides = int(v)
            elif v is True:
                # unlabeled option
                values = k.split(",")
                if len(values) == 2:
                    sides, rows = [int(x) for x in values]
                elif len(values) == 1:
                    sides = int(values[0])
                else:
                    sys.exit(1)

        logger.info(
            "Prism ice with {0} sides and {1} rows.".format(sides, rows))

        L = 2.75
        self.bondlen = 3
        R = L / 2 / sin(pi / sides)
        self.density = sides * rows / \
            (L**3 * 400 * rows) * 18 / 6.022e23 * 1e24

        self.waters = []
        self.fixed = []
        for j in range(rows):
            for i in range(sides):
                x = R * cos(i * pi * 2 / sides)
                y = R * sin(i * pi * 2 / sides)
                z = j * L
                self.waters.append([x, y, z])
            for i in range(sides):
                p = j*sides + i
                q = j*sides + (i+1) % sides
                if j%2 == 0:
                    self.fixed.append([p,q])
                else:
                    self.fixed.append([q,p])
                q = ((j+1) % rows)*sides + i
                if i%2 == 0:
                    self.fixed.append([p,q])
                else:
                    self.fixed.append([q,p])
        self.pairs = self.fixed
        logger.debug(self.fixed)

        self.coord = "absolute"

        self.cell = cellvectors(a=L * 20,
                                b=L * 20,
                                c=L * rows)
