# coding: utf-8
"""
Hydrogen bond network in @NGPH format.
"""

import numpy as np


def hook2(lattice):
    lattice.logger.info("Hook2: Output the undirected network.")

    s = ""
    s += "@NGPH\n"
    s += "{0}\n".format(len(lattice.reppositions))
    for i,j,k in lattice.graph.edges(data=True):
        s += "{0} {1}\n".format(i,j)
    s += "-1 -1\n"
    s = "\n".join(lattice.doc) + "\n" + s
    print(s,end="")
    lattice.logger.info("Hook2: end.")

hooks = {2:hook2}
