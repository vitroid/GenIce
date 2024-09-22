# coding: utf-8

from genice2.cell import cellvectors
import genice2.lattices
desc = {"ref": {},
        "usage": "No options available.",
        "brief": "Conventional high-pressure ice VII."
        }


class Lattice(genice2.lattices.Lattice):
    def __init__(self):
        self.density = 1.6  # default self.density

        self.cell = """
        4 4 4
        """

        self.waters = """
        0.00 0.00 0.00
        0.50 0.50 0.00
        0.50 0.00 0.50
        0.00 0.50 0.50
        0.25 0.25 0.25
        0.75 0.75 0.25
        0.75 0.25 0.75
        0.25 0.75 0.75
        0.50 0.00 0.00
        0.00 0.50 0.00
        0.00 0.00 0.50
        0.50 0.50 0.50
        0.75 0.25 0.25
        0.25 0.75 0.25
        0.25 0.25 0.75
        0.75 0.75 0.75
        """
        self.coord = "relative"
        double_network = True  # It is necessary only for ices 6 and 7

        self.pairs = """
        0 4
        0 6
        1 4
        1 5
        1 6
        2 6
        3 4
        3 5
        3 6
        3 7
        4 2
        5 0
        5 2
        7 0
        7 1
        7 2
        8 12
        8 14
        9 12
        9 13
        9 14
        10 14
        11 12
        11 13
        11 14
        11 15
        12 10
        13 8
        13 10
        15 8
        15 9
        15 10
        """

        self.cell = cellvectors(a=4.0,
                                b=4.0,
                                c=4.0)
