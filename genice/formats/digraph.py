# coding: utf-8
"""
Hydrogen bond network in @NGPH format.
"""

import numpy as np


def hook4(lattice):
    lattice.logger.info("Hook4: Output the hydrogen bond network.")

    s = ""
    s += "@NGPH\n"
    s += "{0}\n".format(len(lattice.reppositions))
    for i,j,k in lattice.spacegraph.edges(data=True):
        s += "{0} {1}\n".format(i,j)
    s += "-1 -1\n"
    s = "\n".join(lattice.doc) + "\n" + s
    print(s,end="")
    lattice.logger.info("Hook4: end.")

hooks = {4:hook4}
