# coding: utf-8

desc={"ref": {"H2": "https://www.britannica.com/science/hydrogen"},
      "brief": "Hydrogen molecule.",
      "usage": "No options available."
      }

import numpy as np
from logging import getLogger
import genice2.molecules

class Molecule(genice2.molecules.Molecule):
    def __init__(self):
        self.sites = np.array([[0,0,-0.037],
                               [0,0,+0.037]])   # nm, HH

        self.labels = ["H","H"]
        self.name = "H2"
