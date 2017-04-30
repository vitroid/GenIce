import logging

import numpy as np

def run(lattice, water_type="TIP3P", guests=[]):
    """
    Hydrogen bond network in @NGPH format.
    """
    logger = logging.getLogger()
    lattice.stage1()   #replicate the unit cell
    lattice.stage2()   #prepare random graph
    logger.info("Output the proximity network.")

    s = ""
    s += "@NGPH\n"
    s += "{0}\n".format(len(lattice.reppositions))
    for i,j,k in lattice.graph.edges_iter(data=True):
        s += "{0} {1}\n".format(i,j)
    s += "-1 -1\n"
    s = "\n".join(lattice.doc) + "\n" + s
    print(s,end="")
