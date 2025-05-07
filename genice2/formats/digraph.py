# coding: utf-8

from genice2.decorators import timeit, banner
import genice2.formats
from logging import getLogger
import sys
from io import TextIOWrapper
from genice2.genice import GenIce

desc = {
    "ref": {"NGPH": "https://vitroid.github.io/@NGPH"},
    "brief": "Directed graph of HBs.",
    "usage": "No options available.",
}


class Format(genice2.formats.Format):
    """
    The topology of the hydrogen bond network is output as a digraph in @NGPH format.
    No options available.
    """

    def __init__(self, **kwargs):
        pass

    @timeit
    @banner
    def dump(self, genice: GenIce, file: TextIOWrapper):
        "Output the hydrogen bond network."
        logger = getLogger()
        s = ""
        s += "@NGPH\n"
        s += "{0}\n".format(len(genice.water_positions()))
        for i, j, k in genice.hydrogen_bond_digraph().edges(data=True):
            s += "{0} {1}\n".format(i, j)
        s += "-1 -1\n"
        s = "\nCommand line: " + " ".join(sys.argv) + "\n" + s
        file.write(s)
