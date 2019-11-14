# coding: utf-8
desc={"ref": {},
      "usage": 'genice eleven[hcchchcc]; Specify layer types with "c" or "h".',
      "brief": "Ice XI w/ stacking faults."
      }

import logging
from genice.cell import cellvectors
import numpy as np

    
lat = [[[0,0], [2,0], [1,3], [3,3]],
       [[0,2], [2,2], [1,5], [3,5]],
       [[0,4], [2,4], [1,1], [3,1]]]

def argparser(arg):
    logger = logging.getLogger()
    global bondlen, density, cell, waters, coord, fixed, pairs
    layer = 0
    height = 0
    dir = 1
    L = []
    N = 0   # number of molecules
    edges = [] # hydrogen bonds to be fixed
    for ch in arg:
        grid = dict()
        right, left = None, None
        for x,y in lat[layer]:
            L.append([x, y, height])
            if right is None:
                right = (x,y)
            grid[x,y] = N
            N += 1
        layer = (layer+dir+3) % 3
        height += 1
        for x,y in lat[layer]:
            L.append([x, y, height])
            if left is None:
                left = (x,y)
            grid[x,y] = N
            N += 1
        height += 3
        # connection along x axis
        x,y = right
        dy = +1
        if (x+1,y+dy) not in grid:
            dy = -1
        for i in range(4):
            A = grid[x,y]
            y0 = (y-dy*2+6) % 6
            C = grid[x,y0]
            edges.append((A,C))
            x = (x+1) % 4
            y = (y+dy+6) % 6
            B = grid[x,y]
            edges.append((A,B))
            dy = -dy
        x,y = left
        dy = -dy
        for i in range(4):
            A = grid[x,y]
            x = (x-1+4) % 4
            y = (y+dy+6) % 6
            B = grid[x,y]
            edges.append((A,B))
            dy = -dy
        # connection along z
        edges.append((N-4, N))
        edges.append((N-3, N+1))
        edges.append((N+2, N-2))
        edges.append((N+3, N-1))
        assert ch in "CcHh"
        if ch in "Hh":
            # hexagonal = alternative
            dir = -dir
            #cubic = progressive
    assert layer == 0 and dir == 1, "Incompatible number of layers."
    fixed = [(A%N, B%N) for A,B in edges]
    pairs = fixed.copy()
    waters = np.array(L) / np.array([4.0, 6.0, height])
    coord = "relative"
    logger.info(waters.shape)
    LHB = 0.276
    bondlen = 0.3
    y = LHB* (8**0.5 / 3)*3
    x = y * 2 / 3**0.5
    z = LHB*height/3
    cell = cellvectors(x,y,z)
    density=0.92

# default.
argparser("hh")
