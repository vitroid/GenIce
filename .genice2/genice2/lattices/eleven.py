# coding: utf-8
"""
Generate a hydrogen-ordered ice XI with stacking faults.

Usage:
  genice2 eleven[hcchchcc]            Specify layer types with "c" or "h".

"""


import genice2.lattices
import numpy as np
from genice2.cell import cellvectors
from logging import getLogger


def usage():
    logger = getLogger()
    logger.info(__doc__)


desc = {
    "ref": {},
    "usage": usage(),
    "brief": "Ice XI w/ stacking faults.",
    "test": ("[ccchchc]",)
}


lat = [[[0, 0], [2, 0], [1, 3], [3, 3]],
       [[0, 2], [2, 2], [1, 5], [3, 5]],
       [[0, 4], [2, 4], [1, 1], [3, 1]]]


class Lattice(genice2.lattices.Lattice):
    def __init__(self, **kwargs):
        logger = getLogger()
        arg = 'hh'
        for k, v in kwargs.items():
            if k == 'layers':
                arg = v
            elif v:  # in case only the char string is given
                arg = k
            else:
                logger.error(f"Unknown option for one plugin: {k}={v}")
        layer = 0
        height = 0
        dir = 1
        L = []
        N = 0   # number of molecules
        edges = []  # hydrogen bonds to be fixed
        for ch in arg:
            grid = dict()
            right, left = None, None
            for x, y in lat[layer]:
                L.append([x, y, height])
                if right is None:
                    right = (x, y)
                grid[x, y] = N
                N += 1
            layer = (layer + dir + 3) % 3
            height += 1
            for x, y in lat[layer]:
                L.append([x, y, height])
                if left is None:
                    left = (x, y)
                grid[x, y] = N
                N += 1
            height += 3
            # connection along x axis
            x, y = right
            dy = +1
            if (x + 1, y + dy) not in grid:
                dy = -1
            for i in range(4):
                A = grid[x, y]
                y0 = (y - dy * 2 + 6) % 6
                C = grid[x, y0]
                edges.append((A, C))
                x = (x + 1) % 4
                y = (y + dy + 6) % 6
                B = grid[x, y]
                edges.append((A, B))
                dy = -dy
            x, y = left
            dy = -dy
            for i in range(4):
                A = grid[x, y]
                x = (x - 1 + 4) % 4
                y = (y + dy + 6) % 6
                B = grid[x, y]
                edges.append((A, B))
                dy = -dy
            # connection along z
            edges.append((N - 4, N))
            edges.append((N - 3, N + 1))
            edges.append((N + 2, N - 2))
            edges.append((N + 3, N - 1))
            assert ch in "CcHh"
            if ch in "Hh":
                # hexagonal = alternative
                dir = -dir
                #cubic = progressive
        assert layer == 0 and dir == 1, "Incompatible number of layers."
        self.fixed = [(A % N, B % N) for A, B in edges]
        self.pairs = self.fixed.copy()
        self.waters = np.array(L) / np.array([4.0, 6.0, height])
        self.coord = "relative"
        LHB = 0.276
        self.bondlen = 0.3
        y = LHB * (8**0.5 / 3) * 3
        x = y * 2 / 3**0.5
        z = LHB * height / 3
        self.cell = cellvectors(x, y, z)
        self.density = 0.92
