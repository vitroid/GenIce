# coding: utf-8
"""
towhee file format (.coods?)
"""

import numpy as np

from genice import rigid

def hook7(lattice):
    lattice.logger.info("Hook7: Output in TowHee format.")
    cellmat = lattice.repcell.mat
    s = ""
    for i,atom in enumerate(lattice.atoms):
        molorder, resname, atomname, position, order = atom
        s += "{0:20.15f}{1:20.15f}{2:20.15f}   {3:6}{4:3}{5:6}{6:6}\n".format(position[0]*10,position[1]*10,position[2]*10,atomname,resname,order+1,i+1)
    print(s,end="")
    s = "  Cell Dimensions:"
    lattice.logger.info(s)
    s = "    a: {0} {1} {2}".format(cellmat[0,0]*10,
                                       cellmat[0,1]*10,
                                       cellmat[0,2]*10)
    lattice.logger.info(s)    
    s = "    b: {0} {1} {2}".format(cellmat[1,0]*10,
                                       cellmat[1,1]*10,
                                       cellmat[1,2]*10)
    lattice.logger.info(s)
    s = "    c: {0} {1} {2}".format(cellmat[2,0]*10,
                                       cellmat[2,1]*10,
                                       cellmat[2,2]*10)
    lattice.logger.info(s)
    lattice.logger.info("Hook7: end.")


hooks = {7:hook7}
