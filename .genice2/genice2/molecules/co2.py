# coding: utf-8
import numpy as np

from logging import getLogger
import genice2.molecules


class Molecule(genice2.molecules.Molecule):
    def __init__(self):
        self.sites_ = np.array(
            [[0, 0, 0], [0, 0, -0.1149], [0, 0, +0.1149]]
        )  # nm, OHHM

        self.labels_ = ["C", "O", "O"]
        self.name_ = "CO2"
