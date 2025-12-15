# coding: utf-8
import math
import numpy as np

from logging import getLogger
import genice2.molecules

water = 1  # Identify


class Molecule(genice2.molecules.Molecule):

    def __init__(self):
        oh = 0.10000
        hangle = 109.47 * math.pi / 180 / 2
        mass = 18
        ohz = oh * math.cos(hangle)
        ohy = oh * math.sin(hangle)
        oz = -ohz * 2 / mass

        self.sites_ = np.array(
            [[0, 0, oz], [0, ohy, ohz + oz], [0, -ohy, ohz + oz]]
        )  # nm

        self.labels_ = ["Ow", "Hw", "Hw"]
        self.name_ = "SOL"
