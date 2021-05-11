# coding: utf-8

from math import pi, sin, cos
import genice2.molecules
from logging import getLogger
import numpy as np
desc = {
    "ref": {"TIP4P(a)": "Jorgensen 1983",
            "TIP4P(b)": "Jorgensen 1985"},
    "usage": "No options available.",
    "brief": "A typical 4-site model.",
}
water = 1  # Identify


class Molecule(genice2.molecules.Molecule):
    def __init__(self):
        L1 = 0.9572 / 10
        L2 = 0.15 / 10
        theta = 104.52 * pi / 180

        hy = L1 * sin(theta / 2)
        hz = L1 * cos(theta / 2)
        mz = L2
        self.sites_ = np.array([[0.0, 0.0, 0.0],
                               [0.0, hy, hz],
                               [0.0, -hy, hz],
                               [0.0, 0.0, mz]])
        self.sites_ -= (self.sites_[1] +
                        self.sites_[2] + self.sites_[3] * 0) / 18
        self.labels_ = ["OW", "HW1", "HW2", "MW"]
        self.name_ = "ICE"
