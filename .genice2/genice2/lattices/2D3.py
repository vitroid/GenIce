from genice2.cell import cellvectors
import genice2.lattices
import numpy as np
desc = {"ref": {},
        "usage": "No options available.",
        "brief": "Trilayer honeycomb ice."
        }


# Unnecessary when pairs are given.
#bondlen = 2.0 * 3.0 / 8.0**0.5 * 1.1


class Lattice(genice2.lattices.Lattice):
    def __init__(self):
        self.cell = [4 * 3**0.5, 6.0, 2 / 8**0.5 * 20]

        # Final self.density of ice.
        self.density = 0.75

        self.coord = "relative"

        self.waters = """
        0.000	0.000	0.
        0.750	0.167	0.05
        0.250	0.167	0.05
        0.750	0.500	0.
        0.250	0.500	0.
        0.500	0.667	0.05
        0.000	0.667	0.05
        0.500	0.000	0.
        0.000	0.000	0.25
        0.750	0.167	0.20
        0.250	0.167	0.20
        0.750	0.500	0.25
        0.250	0.500	0.25
        0.500	0.667	0.20
        0.000	0.667	0.20
        0.500	0.000	0.25
        0.000	0.000	0.40
        0.750	0.167	0.45
        0.250	0.167	0.45
        0.750	0.500	0.40
        0.250	0.500	0.40
        0.500	0.667	0.45
        0.000	0.667	0.45
        0.500	0.000	0.40
        """

        self.pairs = """
        1 0
        1 3
        1 7
        3 5
        3 6
        4 6
        4 5
        5 7
        6 0
        0 2
        7 2
        4 2
        1 9
        5 13
        6 14
        2 10
        9 8
        9 11
        9 15
        11 13
        11 14
        12 14
        12 13
        13 15
        14 8
        8 10
        15 10
        12 10
        8 16
        11 19
        12 20
        15 23
        17 16
        17 19
        17 23
        19 21
        19 22
        20 22
        20 21
        21 23
        22 16
        16 18
        23 18
        20 18
        """

        self.cell = cellvectors(a=6.928203230275509,
                                b=6.0,
                                c=14.14213562373095)
