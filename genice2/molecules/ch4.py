# coding: utf-8
import numpy as np

from logging import getLogger
import genice2.molecules


desc = {"usage": "No options available.", "brief": "An all-atom methane model."}


class Molecule(genice2.molecules.Molecule):
    def __init__(self):
        self.sites_ = np.array(
            [
                [0.0, 0.0, 0.0],
                [-1.0, -1.0, -1.0],
                [-1.0, +1.0, +1.0],
                [+1.0, -1.0, +1.0],
                [+1.0, +1.0, -1.0],
            ]
        )  # CHHHH
        CH = 0.109  # nm
        self.sites_ *= CH / (3.0**0.5)

        self.labels_ = ["C", "H", "H", "H", "H"]
        self.name_ = "CH4"
