# coding: utf-8
"""
MDView file format
"""

import numpy as np


def hook7(lattice):
    lattice.logger.info("Hook7: Output in MDView format.")
    lattice.logger.info("  Total number of atoms: {0}".format(len(lattice.atoms)))
    s = ""
    #if celltype == "rect":
    #    s += "-length '({0}, {1}, {2})'\n".format(repcell[0,0]*10,repcell[1,1]*10,repcell[2,2]*10)
    s += "-center 0 0 0\n"
    s += "-fold\n"
    s += "{0}\n".format(len(lattice.atoms))
    for i in range(len(lattice.atoms)):
        molorder, resname, atomname, position, order = lattice.atoms[i]
        s += "{0:5} {1:9.4f} {2:9.4f} {3:9.4f}\n".format(atomname,position[0]*10,position[1]*10,position[2]*10)
    s = '#' + "\n#".join(lattice.doc) + "\n" + s
    print(s,end="")
    lattice.logger.info("Hook7: end.")


hooks = {7:hook7}
