# coding: utf-8

from genice2.decorators import timeit, banner
import genice2.formats
from logging import getLogger
desc = {"ref": {"NGPH": "https://vitroid.github.io/@NGPH"},
        "brief": "Directed graph of HBs.",
        "usage": "No options available."
        }


class Format(genice2.formats.Format):
    """
The topology of the hydrogen bond network is output as a digraph in @NGPH format.
No options available.
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def hooks(self):
        return {4: self.Hook4}

    @timeit
    @banner
    def Hook4(self, ice):
        "Output the hydrogen bond network."
        logger = getLogger()
        s = ""
        s += "@NGPH\n"
        s += "{0}\n".format(len(ice.reppositions))
        for i, j, k in ice.spacegraph.edges(data=True):
            s += "{0} {1}\n".format(i, j)
        s += "-1 -1\n"
        s = "\n".join(ice.doc) + "\n" + s
        self.output = s
