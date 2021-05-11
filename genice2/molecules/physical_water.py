# coding: utf-8

import genice2.molecules
from logging import getLogger
import numpy as np
import math
desc = {
    "ref": {
        "TIP3P": "Jorgensen, Chandrasekhar, Madura, Impey, Klein, J Chem Phys, 79, 926 (1983)."},
    "brief": "Physical model of water; Oxygen atom is on the lattice point.",
    "usage": "No options available."}


water = 1  # Identify


class Molecule(genice2.molecules.Molecule):
    def __init__(self):
        oh = 0.09572
        hangle = 104.52 * math.pi / 180 / 2
        mass = 18
        ohz = oh * math.cos(hangle)
        ohy = oh * math.sin(hangle)
        oz = 0  # -ohz*2/mass
        self.sites_ = np.array([[0, 0, oz],
                                [0, ohy, ohz + oz],
                                [0, -ohy, ohz + oz]])  # nm, OHHM

        self.labels_ = ["O", "H", "H"]
        self.name_ = "SOL"
