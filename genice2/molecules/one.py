# coding: utf-8
import numpy as np
import genice2.molecules
from dataclasses import dataclass


# @dataclass
class Molecule(genice2.molecules.Molecule):
    def __init__(self, label: str = "X", name: str = "X"):
        self.sites = np.array([[0, 0, 0]])
        self.labels = [
            label,
        ]
        self.name = name
