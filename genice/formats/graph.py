# coding: utf-8
desc={"ref": {"NGPH": "https://vitroid.github.io/@NGPH"},
      "brief": "Undirected graph of HBs.",
      "usage": "No options available."
      }



import numpy as np
from logging import getLogger

import genice.formats
class Format(genice.formats.Format):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)


    def hooks(self):
        return {2:self.hook2}


    def hook2(self, lattice):
        logger = getLogger()
        logger.info("Hook2: Output the undirected network.")

        s = ""
        s += "@NGPH\n"
        s += "{0}\n".format(len(lattice.reppositions))
        for i,j,k in lattice.graph.edges(data=True):
            s += "{0} {1}\n".format(i,j)
        s += "-1 -1\n"
        s = "\n".join(lattice.doc) + "\n" + s
        print(s,end="")
        logger.info("Hook2: end.")
