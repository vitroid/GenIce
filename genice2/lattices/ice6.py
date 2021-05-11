# coding: utf-8

from genice2.cell import cellvectors
import genice2.lattices
desc = {"ref": {"VI": 'Petrenko 1999'},
        "usage": "No options available.",
        "brief": "Conventional high-pressure ice VI."
        }


class Lattice(genice2.lattices.Lattice):
    def __init__(self):
        self.density = 1.373  # default self.density

        self.bondlen = 3  # bond threshold
        self.cell = """
        6.181 6.181 5.698
        """

        double_network = True  # It is necessary only for ices 6 and 7

        self.coord = "relative"

        self.waters = """
            0.7504    0.2502    0.7510
            0.4707    0.2502    0.3666
            0.9710    0.7509    0.6348
            0.5298    0.7509    0.6348
            0.2501    0.4710    0.8674
            0.0295    0.2502    0.3666
            0.7504    0.5301    0.1341
            0.2501    0.7509    0.2503
            0.7504    0.9716    0.1341
            0.2501    0.0295    0.8674
        """

        self.pairs = """
        1 0
        3 2
        4 2
        4 3
        5 0
        5 1
        6 0
        6 1
        6 5
        7 2
        7 3
        7 4
        8 0
        8 1
        8 5
        8 6
        9 2
        9 3
        9 4
        9 7
        """

        self.cell = cellvectors(a=6.181,
                                b=6.181,
                                c=5.698)
