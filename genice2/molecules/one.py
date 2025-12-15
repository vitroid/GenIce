# coding: utf-8
import numpy as np
from logging import getLogger
import genice2.molecules


class Molecule(genice2.molecules.Molecule):
    def __init__(self, **kwargs):
        # placeholder for 12-hedral cage
        self.sites_ = np.array([[0, 0, 0]])
        self.labels_ = [
            "X",
        ]
        self.name_ = "X"
        if "label" in kwargs:
            self.labels_ = [
                kwargs["label"],
            ]
        if "name" in kwargs:
            self.name_ = kwargs["name"]
