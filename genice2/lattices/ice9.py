# coding: utf-8

from genice2.cell import cellvectors
import genice2.lattices
desc = {"ref": {"9": 'Londono 1993'},
        "usage": "No options available.",
        "brief": "Ice IX, a hydrogen-ordered counterpart of ice III."
        }


class Lattice(genice2.lattices.Lattice):
    def __init__(self):
        self.waters = """
        0.895 0.105 0.750
        0.105 0.895 0.250
        0.605 0.605 0.500
        0.395 0.395 0.000
        0.298 0.108 0.710
        0.702 0.892 0.210
        0.392 0.798 0.960
        0.608 0.202 0.460
        0.798 0.392 0.040
        0.202 0.608 0.540
        0.892 0.702 0.790
        0.108 0.298 0.290
        """
        self.coord = "relative"
        self.bondlen = 3
        self.density = 1.15672

        self.fixed = """
        0 7
        0 8
        1 6
        1 9
        2 5
        2 10
        3 4
        3 11
        4 0
        4 6
        5 1
        5 7
        6 5
        6 3
        7 2
        7 4
        8 10
        8 3
        10 0
        10 9
        9 2
        9 11
        11 1
        11 8
        """

        self.cell = cellvectors(a=6.67,
                                b=6.67,
                                c=6.97)
