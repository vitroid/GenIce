# coding: utf-8
desc={"ref": {"TIP4P": "W. L. Jorgensen, J. Chandrasekhar, J. D. Madura, R. W. Impey, and M. L. Klein, "
+"Comparison of simple potential functions for simulating liquid water, J. Chem. Phys. 79 (1983) 926-935; "
+"W. L. Jorgensen and J. D. Madura, Temperature and size dependence for monte carlo simulations of TIP4P water, "
+"Mol. Phys. 56 (1985) 1381-1392."},
      "brief": "A typical four-site water model.",
      "usage": "No options available."
      }


import numpy as np
from math import pi,sin,cos

L1 = 0.9572 / 10
L2 = 0.15 / 10
theta=104.52 * pi/180


hy = L1*sin(theta/2)
hz = L1*cos(theta/2)
mz = L2

sites = np.array([[0.0, 0.0, 0.0],
                  [0.0, hy,  hz],
                  [0.0,-hy,  hz],
                  [0.0, 0.0, mz]])
sites -= (sites[1]+sites[2]+sites[3]*16)/18


atoms = ["O","H","H","."]
labels = ["OW","HW1","HW2","MW"]
name = "ICE"
