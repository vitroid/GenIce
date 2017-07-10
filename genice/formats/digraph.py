# coding: utf-8
import logging

import numpy as np
from genice import formatter

class Formatter(formatter.Formatter):
    """
    Hydrogen bond network in @NGPH format.
    """
    def __init__(self):
        self.hooks[4] = self.hook4

        
    def hook4(self, lattice):
        logger = logging.getLogger()
        logger.info("Output the hydrogen bond network.")
    
        s = ""
        s += "@NGPH\n"
        s += "{0}\n".format(len(lattice.reppositions))
        for i,j,k in lattice.graph.edges_iter(data=True):
            s += "{0} {1}\n".format(i,j)
        s += "-1 -1\n"
        s = "\n".join(lattice.doc) + "\n" + s
        print(s,end="")
