# coding: utf-8
import numpy as np

# All-atom THF model

from logging import getLogger
import genice2.molecules


desc = {
    "usage": "No options available.",
    "brief": "An all-atom tetrahydrofuran (THF) model.",
}


class Molecule(genice2.molecules.Molecule):
    def __init__(self):
        self.sites_ = (
            np.array(
                [
                    [1.2328, -0.0005, 0.0000],
                    [-1.0107, -0.7202, -0.2205],
                    [-1.0102, 0.7210, 0.2205],
                    [0.3936, -1.1560, 0.1374],
                    [0.3946, 1.1557, -0.1375],
                    [-1.7823, -1.3279, 0.2593],
                    [-1.1544, -0.7757, -1.3060],
                    [-1.7812, 1.3292, -0.2593],
                    [-1.1537, 0.7766, 1.3061],
                    [0.4518, -1.4889, 1.1792],
                    [0.7622, -1.9589, -0.5071],
                    [0.4532, 1.4885, -1.1793],
                    [0.7639, 1.9583, 0.5070],
                ]
            )
            / 10
        )

        self.atoms_ = ["O"] + ["C"] * 4 + ["H"] * 8
        self.labels_ = ["O", "CA", "CA", "CB", "CB"] + ["H"] * 8
        self.name_ = "THF"
