# coding: utf-8

from genice2.cell import cellvectors
import genice2.lattices
desc = {
    "ref": {
        "sVII": "Jeffrey 1984",
        "CS4": "Kosyakov 1999",
        "SOD": "IZA Database",
        "207_1_4435": "Engel 2018",
        "engel01": "Engel 2018",
    },
    "usage": "No options available.",
    "brief": "Hypothetical zeolitic ice"
}


class Lattice(genice2.lattices.Lattice):
    def __init__(self):
        self.cell = """
        7.8064588643 7.8064588643 7.8064588643
        """
        self.waters = """
        0.0 5.85484414822 3.90322943215
        3.90322943215 1.95161471607 0.0
        1.95161471607 0.0 3.90322943215
        5.85484414822 3.90322943215 0.0
        0.0 3.90322943215 1.95161471607
        3.90322943215 0.0 5.85484414822
        0.0 3.90322943215 5.85484414822
        3.90322943215 0.0 1.95161471607
        3.90322943215 5.85484414822 0.0
        0.0 1.95161471607 3.90322943215
        5.85484414822 0.0 3.90322943215
        1.95161471607 3.90322943215 0.0
        """
        self.coord = "absolute"
        self.bondlen = 3
        self.density = 0.75396428378

        self.cages = """
        K 0 0 0
        K 0.5 0.5 0.5
        """

        self.cell = cellvectors(a=7.8064588643,
                                b=7.8064588643,
                                c=7.8064588643)
