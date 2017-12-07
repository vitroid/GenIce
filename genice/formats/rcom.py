# coding: utf-8
"""
Centers of mass of water molecule 
"""

import numpy as np


def hook1(lattice):
    lattice.logger.info("Hook1: Output centers of mass of water molecules.")
    cellmat = lattice.repcell.mat
    s = ""
    if cellmat[1,0] == 0 and cellmat[2,0] == 0 and cellmat[2,1] == 0:
        s += "@BOX3\n"
        s += "{0} {1} {2}\n".format(cellmat[0,0]*10,cellmat[1,1]*10,cellmat[2,2]*10)
    else:
        s += "@BOX9\n"
        for d in range(3):
            s += "{0} {1} {2}\n".format(cellmat[0,d]*10,cellmat[1,d]*10,cellmat[2,d]*10)
    s += "@AR3R\n"
    s += "{0}\n".format(len(lattice.reppositions))
    for rpos in lattice.reppositions:
        s += "{0:9.4f} {1:9.4f} {2:9.4f}\n".format(rpos[0],
                                                rpos[1],
                                                rpos[2])
    s = "\n".join(lattice.doc) + "\n" + s
    print(s,end="")
    lattice.logger.info("Hook1: end.")


hooks = {1:hook1}
