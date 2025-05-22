# coding: utf-8
"""
Centers of mass of water molecule
"""

from logging import getLogger
from io import TextIOWrapper
import sys

import numpy as np

from genice2.decorators import timeit, banner
import genice2.formats
from genice2.genice import GenIce


class Format(genice2.formats.Format):
    def __init__(self, **kwargs):
        pass

    @timeit
    @banner
    def dump(self, genice: GenIce, file: TextIOWrapper):
        "Output centers of mass of water molecules."
        logger = getLogger()
        cellmat = genice.cell_matrix()
        s = ""
        if cellmat[1, 0] == 0 and cellmat[2, 0] == 0 and cellmat[2, 1] == 0:
            s += "@BOX3\n"
            s += "{0} {1} {2}\n".format(
                cellmat[0, 0] * 10, cellmat[1, 1] * 10, cellmat[2, 2] * 10
            )
        else:
            s += "@BOX9\n"
            for d in range(3):
                s += "{0} {1} {2}\n".format(
                    cellmat[0, d] * 10, cellmat[1, d] * 10, cellmat[2, d] * 10
                )
        s += "@AR3R\n"
        water_positions = genice.water_positions()
        s += "{0}\n".format(len(water_positions))
        for rpos in water_positions:
            s += "{0:9.4f} {1:9.4f} {2:9.4f}\n".format(rpos[0], rpos[1], rpos[2])
        s = "\nCommand line: " + " ".join(sys.argv) + "\n" + s
        file.write(s)
