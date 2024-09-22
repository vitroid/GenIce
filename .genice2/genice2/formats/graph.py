# coding: utf-8
import genice2.formats
from genice2.decorators import timeit, banner
from logging import getLogger
import numpy as np
desc = {"ref": {"NGPH": "https://vitroid.github.io/@NGPH"},
        "brief": "Undirected graph of HBs.",
        "usage": "No options available."
        }


class Format(genice2.formats.Format):
    """
The topology of the hydrogen bond network is output as an undirected graph in @NGPH format.
No options available.
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def hooks(self):
        return {2: self.Hook2}

    @timeit
    @banner
    def Hook2(self, ice):
        "Output the undirected network."
        logger = getLogger()

        s = ""
        s += "@NGPH\n"
        s += "{0}\n".format(len(ice.reppositions))
        for i, j, k in ice.graph.edges(data=True):
            s += "{0} {1}\n".format(i, j)
        s += "-1 -1\n"
        s = "\n".join(ice.doc) + "\n" + s
        self.output = s
