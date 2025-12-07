# coding: utf-8
import numpy as np
import genice3.molecule


class Molecule(genice3.molecule.Molecule):
    def __init__(self, label: str = "X", name: str = "X"):
        self.sites = np.array([[0, 0, 0]])
        self.labels = [
            label,
        ]
        self.name = name
