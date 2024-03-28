# coding: utf-8
import numpy as np

# United-atom THF model

from logging import getLogger
import genice2.molecules


desc = {
    "usage": "No options available.",
    "brief": "A united-atom five-site tetrahydrofuran (THF) model.",
}


class Molecule(genice2.molecules.Molecule):
    def __init__(self):
        self.sites_ = np.array(
            [
                [0.0, -0.119625, 0.0],
                [0.116284, -0.039705, 0.0],
                [0.076453, 0.107915, 0.0],
                [-0.076453, 0.107915, 0.0],
                [-0.116284, -0.039705, 0.0],
            ]
        )

        self.atoms_ = ["O", "C", "C", "C", "C"]
        self.labels_ = ["O", "CA", "CB", "CB", "CA"]
        self.name_ = "THF"
