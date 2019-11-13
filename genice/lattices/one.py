# coding: utf-8
desc={"ref": {},
      "usage": "one[hcchchcc] last layer type is determined automatically.",
      "brief": "Ice I w/ stacking faults."
      }

import logging
from genice.cell import cellvectors
import numpy as np

    
layers = np.array([[[1/4, 0, 0], [3/4, 0, 0], [0, 1/2, 0], [1/2, 1/2, 0]],
                   [[0, 5/6, 0], [1/2, 5/6, 0], [1/4, 1/3, 0], [3/4, 1/3, 0]],
                   [[1/4, 2/3, 0], [3/4, 2/3, 0], [0,1/6, 0], [1/2, 1/6, 0]]])

def argparser(arg):
    logger = logging.getLogger()
    global sides, rows, bondlen, density, cell, waters, coord
    layer = 0
    height = 0
    dir = 1
    L = []
    L.append(layers[layer] + np.array([0.0, 0.0, height]))
    layer = (layer+dir+3) % 3
    height += 1
    L.append(layers[layer] + np.array([0.0, 0.0, height]))
    height += 3
    for ch in arg:
        if ch in "Hh":
            # hexagonal = alternative
            dir = -dir
        else:
            #cubic = progressive
            pass
        L.append(layers[layer] + np.array([0.0, 0.0, height]))
        layer = (layer+dir+3) % 3
        height += 1
        L.append(layers[layer] + np.array([0.0, 0.0, height]))
        height += 3
    assert layer == 0, "Incompatible number of layers."
    L = L / np.array([1.0, 1.0, height])
    waters = L.reshape(L.shape[0]*4, 3)
    logger.info(waters.shape)
    LHB = 0.276
    bondlen = 0.3
    y = LHB* (8**0.5 / 3)*3
    x = y * 2 / 3**0.5
    z = LHB*height/3
    coord = "relative"
    cell = cellvectors(x,y,z)
    density=0.92

# default.
argparser("hhcc")
