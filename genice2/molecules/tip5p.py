# coding: utf-8
import math
import numpy as np
from logging import getLogger
import genice2.molecules

desc = {
    "usage": "No options available.",
    "brief": "A typical 5-site model.",
}

water = 1  # Identify


class Molecule(genice2.molecules.Molecule):

    def __init__(self):
        oh = 0.09572  # nm
        om = 0.07  # nm
        hangle = 104.52 * math.pi / 180 / 2
        mangle = 109.47 * math.pi / 180 / 2
        mass = 18
        ohz = oh * math.cos(hangle)
        ohy = oh * math.sin(hangle)
        omz = -om * math.cos(mangle)
        omx = om * math.sin(mangle)
        oz = -ohz * 2 / mass
        self.sites_ = np.array([[0, 0, oz],
                                [0, ohy, ohz + oz],
                                [0, -ohy, ohz + oz],
                                [omx, 0, omz + oz],
                                [-omx, 0, omz + oz],
                                ])  # nm, OHHMM

        self.atoms_ = ["O", "H", "H", ".", "."]
        self.labels_ = ["OW", "HW1", "HW2", "MW1", "MW2"]
        self.name_ = "SOL"
