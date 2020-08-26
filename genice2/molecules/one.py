# coding: utf-8
import numpy as np
from logging import getLogger
import genice2.molecules

class Molecule(genice2.molecules.Molecule):
    def __init__(self, **kwargs):
        #placeholder for 12-hedral cage
        self.sites = np.array([[0,0,0]])
        self.labels=["X",]
        self.name="X"
        if "label" in kwargs:
            self.labels = [kwargs["label"],]
        if "name" in kwargs:
            self.name = kwargs["name"]
