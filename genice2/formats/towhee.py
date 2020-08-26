# coding: utf-8
"""
towhee file format (.coods?) [EXPERIMENTAL]
"""

import numpy as np
from logging import getLogger
from genice2 import rigid


import genice2.formats
class Format(genice2.formats.Format):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)


    def hooks(self):
        return {7:self.hook7}


    def hook7(self, ice):
        logger = getLogger()
        logger.info("Hook7: Output in TowHee format.")
        cellmat = ice.repcell.mat
        s = ""
        for i,atom in enumerate(ice.atoms):
            molorder, resname, atomname, position, order = atom
            s += "{0:20.15f}{1:20.15f}{2:20.15f}   {3:6}{4:3}{5:6}{6:6}\n".format(position[0]*10,position[1]*10,position[2]*10,atomname,resname,order+1,i+1)
        print(s,end="")
        s = "  Cell Dimensions:"
        logger.info(s)
        s = "    a: {0} {1} {2}".format(cellmat[0,0]*10,
                                           cellmat[0,1]*10,
                                           cellmat[0,2]*10)
        logger.info(s)
        s = "    b: {0} {1} {2}".format(cellmat[1,0]*10,
                                           cellmat[1,1]*10,
                                           cellmat[1,2]*10)
        logger.info(s)
        s = "    c: {0} {1} {2}".format(cellmat[2,0]*10,
                                           cellmat[2,1]*10,
                                           cellmat[2,2]*10)
        logger.info(s)
        logger.info("Hook7: end.")
