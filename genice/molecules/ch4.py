# coding: utf-8
import numpy as np

from logging import getLogger
import genice.molecules

class Molecule(genice.molecules.Molecule):
    def __init__(self):
        self.sites = np.array([[0.0, 0.0, 0.0],
                          [-1.,-1.,-1.],
                          [-1.,+1.,+1.],
                          [+1.,-1.,+1.],
                          [+1.,+1.,-1.]]) #CHHHH
        CH = 0.109  # nm
        self.sites *= CH / (3.0**0.5)

        self.labels = ["C","H","H","H","H"]
        self.name = "CH4"
