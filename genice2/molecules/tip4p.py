# coding: utf-8
"""
[1] M. W. Mahoney and W. L. Jorgensen, A five-site model for liquid water and the reproduction of the density anomaly by rigid, nonpolarizable potential functions, J. Chem. Phys. 112 (2000) 8910-8922. [Back]
[2] W. L. Jorgensen, J. Chandrasekhar, J. D. Madura, R. W. Impey, and M. L. Klein, Comparison of simple potential functions for simulating liquid water, J. Chem. Phys. 79 (1983) 926-935
[3] W. L. Jorgensen and J. D. Madura, Temperature and size dependence for monte carlo simulations of TIP4P water, Mol. Phys. 56 (1985) 1381-1392.
"""
import numpy as np
from math import pi,sin,cos
from logging import getLogger
import genice2.molecules

class Molecule(genice2.molecules.Molecule):
    def __init__(self):
        L1 = 0.9572 / 10
        L2 = 0.15 / 10
        theta=104.52 * pi/180

        hy = L1*sin(theta/2)
        hz = L1*cos(theta/2)
        mz = L2
        self.sites = np.array([[0.0, 0.0, 0.0],
                               [0.0, hy,  hz],
                               [0.0,-hy,  hz],
                               [0.0, 0.0, mz]])
        self.sites -= (self.sites[1]+self.sites[2]+self.sites[3]*0)/18
        self.atoms = ["O","H","H","."]
        self.labels = ["OW","HW1","HW2","MW"]
        self.name = "ICE"
