# coding: utf-8
from logging import getLogger
from io import TextIOWrapper
import sys

import numpy as np

import genice2.formats
from genice2.decorators import timeit, banner
from genice2.genice import GenIce

desc = {
    "ref": {"NGPH": "https://vitroid.github.io/@NGPH"},
    "brief": "Undirected graph of HBs.",
    "usage": "No options available.",
}


class Format(genice2.formats.Format):
    """
    The topology of the hydrogen bond network is output as an undirected graph in @NGPH format.
    No options available.
    """

    def __init__(self, **kwargs):
        pass

    @timeit
    @banner
    def dump(self, genice: GenIce, file: TextIOWrapper):
        "Output the undirected network."
        logger = getLogger()

        s = ""
        s += "@NGPH\n"
        s += "{0}\n".format(len(genice.water_positions()))
        for i, j, k in genice.hydrogen_bond_graph().edges(data=True):
            s += "{0} {1}\n".format(i, j)
        s += "-1 -1\n"
        s = "\nCommand line: " + " ".join(sys.argv) + "\n" + s
        file.write(s)
