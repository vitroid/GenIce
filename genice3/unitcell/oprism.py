# coding: utf-8
"""
Generate a hydrogen-ordered prism ice. (Cylindrical ice)

Usage:
    genice2 oprism       A prism of default sides and rows (6,10)
    genice2 oprism[5]    A pentagonal prism of 10 rows.
    genice2 oprism[8,6]  An octagonal prism of 6 rows.

Options:
    [sides[,rows]]    Number of sides and rows.

"""


import genice3.unitcell
from genice3.util import cellvectors
from logging import getLogger
from math import sin, pi, cos
import networkx as nx


desc = {
    "ref": {"prism": "Koga 2001"},
    "usage": __doc__,
    "brief": "Hydrogen-ordered ice nanotubes.",
    "test": [{"args": {}}, {"args": "5"}, {"args": "8,6"}],
}


class UnitCell(genice3.unitcell.UnitCell):
    def __init__(self, **kwargs):
        logger = getLogger()
        # global sides, rows, bondlen, density, cell, waters, coord
        sides = 6
        rows = 10

        for k, v in kwargs.items():
            if k == "rows":
                rows = int(v)
                assert rows % 2 == 0, "Number of rows nums be even."
            elif k == "sides":
                sides = int(v)
            else:
                raise ValueError(f"Unknown option for oprism plugin: {k}={v}")

        logger.info("Prism ice with {0} sides and {1} rows.".format(sides, rows))

        L = 2.75
        bondlen = 3
        R = L / 2 / sin(pi / sides)
        density = sides * rows / (L**3 * 400 * rows) * 18 / 6.022e23 * 1e24

        waters = []
        fixed = []
        for j in range(rows):
            for i in range(sides):
                x = R * cos(i * pi * 2 / sides)
                y = R * sin(i * pi * 2 / sides)
                z = j * L
                waters.append([x, y, z])
            for i in range(sides):
                p = j * sides + i
                q = j * sides + (i + 1) % sides
                if j % 2 == 0:
                    fixed.append([p, q])
                else:
                    fixed.append([q, p])
                q = ((j + 1) % rows) * sides + i
                if i % 2 == 0:
                    fixed.append([p, q])
                else:
                    fixed.append([q, p])

        coord = "absolute"

        cell = cellvectors(a=L * 20, b=L * 20, c=L * rows)

        super().__init__(
            cell=cell,
            waters=waters,
            coord=coord,
            # bondlen=bondlen,
            density=density,
            fixed=nx.DiGraph(fixed),
            graph=nx.Graph(fixed),
        )
