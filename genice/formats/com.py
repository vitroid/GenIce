# coding: utf-8
"""
Centers of mass of water molecule 
"""

import numpy as np

def hook1(lattice):
    lattice.logger.info("Hook1: Output centers of mass of water molecules.")
    s = ""
    if lattice.repcell[1,0] == 0 and lattice.repcell[2,0] == 0 and lattice.repcell[2,1] == 0:
        s += "@BOX3\n"
        s += "{0} {1} {2}\n".format(lattice.repcell[0,0]*10,lattice.repcell[1,1]*10,lattice.repcell[2,2]*10)
    else:
        s += "@BOX9\n"
        for d in range(3):
            s += "{0} {1} {2}\n".format(lattice.repcell[0,d]*10,lattice.repcell[1,d]*10,lattice.repcell[2,d]*10)
    s += "@AR3A\n"
    s += "{0}\n".format(len(lattice.reppositions))
    for i in range(len(lattice.reppositions)):
        position = np.dot(lattice.reppositions[i],lattice.repcell)*10   #in Angstrom
        s += "{0:9.4f} {1:9.4f} {2:9.4f}\n".format(position[0],
                                                position[1],
                                                position[2])
    s = "\n".join(lattice.doc) + "\n" + s
    print(s,end="")
    lattice.logger.info("Hook1: end.")

hooks = {1:hook1}
