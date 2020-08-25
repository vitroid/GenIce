# coding: utf-8
"""
Centers of mass of water molecule
"""

import numpy as np
from logging import getLogger

import genice.formats
class Format(genice.formats.Format):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)


    def hooks(self):
        return {1:self.hook1}


    def hook1(self, ice):
        logger = getLogger()
        logger.info("Hook1: Output centers of mass of water molecules.")
        cellmat = ice.repcell.mat
        s = ""
        if cellmat[1,0] == 0 and cellmat[2,0] == 0 and cellmat[2,1] == 0:
            s += "@BOX3\n"
            s += "{0} {1} {2}\n".format(cellmat[0,0]*10,cellmat[1,1]*10,cellmat[2,2]*10)
        else:
            s += "@BOX9\n"
            for d in range(3):
                s += "{0} {1} {2}\n".format(cellmat[0,d]*10,cellmat[1,d]*10,cellmat[2,d]*10)
        s += "@AR3A\n"
        s += "{0}\n".format(len(ice.reppositions))
        for pos in ice.reppositions:
            position = np.dot(pos,cellmat)*10   #in Angstrom
            s += "{0:9.4f} {1:9.4f} {2:9.4f}\n".format(position[0],
                                                       position[1],
                                                       position[2])
        s = "\n".join(ice.doc) + "\n" + s
        print(s,end="")
        logger.info("Hook1: end.")
