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
    global bondlen, density, cell, waters, coord
    layer = 0
    height = 0
    dir = 1
    L = []
    for ch in arg:
        for x,y in lat[layer]:
            L.append([x, y, height])
        layer = (layer+dir+3) % 3
        height += 1
        for x,y in lat[layer]:
            L.append([x, y, height])
        height += 3
        assert ch in "CcHh"
        if ch in "Hh":
            # hexagonal = alternative
            dir = -dir
            #cubic = progressive
    assert layer == 0 and dir == 1, "Incompatible number of layers."
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
