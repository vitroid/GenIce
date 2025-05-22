# coding: utf-8

import genice2.molecules
from logging import getLogger
import numpy as np

desc = {
    "ref": {"H2": "https://www.britannica.com/science/hydrogen"},
    "brief": "Hydrogen molecule.",
    "usage": "No options available.",
}


class Molecule(genice2.molecules.Molecule):
    def __init__(self):
        self.sites = np.array([[0, 0, -0.037], [0, 0, +0.037]])  # nm, HH
        self.labels = ["H", "H"]
        self.name = "H2"
