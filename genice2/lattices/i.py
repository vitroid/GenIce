# coding: utf-8

import genice2.lattices
from genice2.cell import cellvectors


desc = {"ref": {"0": "Fennell 2005", },
        "usage": "No options available.",
        "brief": 'Hypothetical ice "i".'
        }


class Lattice(genice2.lattices.Lattice):
    def __init__(self):
        self.density = 0.92  # default self.density

        self.bondlen = 0.4  # bond threshold

        self.cell = """
        1.0 1.0 0.94
        """

        self.waters = """
        0.16666 0.16666 0.0
        0.16666 0.16666 0.47
        0.16666 0.83333 0.0
        0.16666 0.83333 0.47
        0.83333 0.16666 0.0
        0.83333 0.16666 0.47
        0.83333 0.83333 0.0
        0.83333 0.83333 0.47
        0.33333 0.33333 0.235
        0.33333 0.33333 0.705
        0.33333 0.66666 0.235
        0.33333 0.66666 0.705
        0.66666 0.33333 0.235
        0.66666 0.33333 0.705
        0.66666 0.66666 0.235
        0.66666 0.66666 0.705
        """

        self.coord = "absolute"

        self.cell = cellvectors(a=1.0,
                                b=1.0,
                                c=0.94)
