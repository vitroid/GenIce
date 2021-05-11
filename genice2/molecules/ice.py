# coding: utf-8
"""
[1] J. L. F. Abascal, E. Sanz, R. G. Fernández, and C. Vega, A potential model for the study of ices and amorphous water: TIP4P/Ice, J. Chem. Phys. 122 (2005) 234511.
"""

from math import pi, sin, cos
import genice2.molecules
from logging import getLogger
import numpy as np
water = 1  # Identify


class Molecule(genice2.molecules.Molecule):

    def __init__(self):
        L1 = 0.9572 / 10
        L2 = 0.1577 / 10
        theta = 104.52 * pi / 180

        hy = L1 * sin(theta / 2)
        hz = L1 * cos(theta / 2)
        mz = L2
        self.sites_ = np.array([[0.0, 0.0, 0.0],
                                [0.0, hy, hz],
                                [0.0, -hy, hz],
                                [0.0, 0.0, mz]])
        self.sites_ -= (self.sites_[1] +
                        self.sites_[2] + self.sites_[0] * 0) / 18

        self.atoms_ = ["O", "H", "H", "."]
        self.labels_ = ["OW", "HW1", "HW2", "MW"]
        self.name_ = "ICE"
