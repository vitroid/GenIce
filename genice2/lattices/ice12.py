# coding: utf-8

from genice2.cell import cellvectors
import genice2.lattices

desc = {
    "ref": {"12(1)": "Lobban 1998", "12(2)": "Koza 2000"},
    "usage": "No options available.",
    "brief": "Metastable high-pressure ice XII.",
}


class Lattice(genice2.lattices.Lattice):
    def __init__(self):
        self.bondlen = 3  # bond threshold

        self.density = 1.4397

        self.waters = """
        0.75 0.13356 0.1875
        0.75 0.13536 0.6875
        0.86464 0.75 0.3125
        0.86464 0.75 0.8125
        0.13536 0.25 0.3125
        0.13536 0.25 0.8125
        0.5 0 0.375
        0.5 0 0.875
        0.25 0.86464 0.1875
        0.25 0.86464 0.6875
        0 0.5 0.125
        0 0.5 0.625
        0.5 0.5 0.25
        0.5 0.5 0.75
        0 0 0
        0 0 0.5
        0.25 0.63536 0.4375
        0.25 0.63536 0.9375
        0.36464 0.25 0.0625
        0.36464 0.25 0.5625
        0.75 0.36564 0.4375
        0.75 0.36564 0.9375
        0.63536 0.75 0.0625
        0.63536 0.75 0.5625
        """

        self.coord = "relative"

        self.cell = cellvectors(a=8.2816, b=8.2816, c=8.0722)
