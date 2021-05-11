# coding: utf-8
import math
import numpy as np
from logging import getLogger
import genice2.molecules

desc = {
    "usage": "No options available.",
    "brief": "A typical 3-site model.",
}
water = 1  # Identify


class Molecule(genice2.molecules.Molecule):
    def __init__(self):
        oh = 0.09572
        hangle = 104.52 * math.pi / 180 / 2
        mass = 18
        ohz = oh * math.cos(hangle)
        ohy = oh * math.sin(hangle)
        oz = -ohz * 2 / mass
        self.sites_ = np.array([[0, 0, oz],
                               [0, ohy, ohz + oz],
                               [0, -ohy, ohz + oz]])  # nm, OHHM
        self.labels_ = ["O", "H", "H"]
        self.name_ = "SOL"
