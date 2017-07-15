# coding: utf-8
"""
[1] J. L. F. Abascal, E. Sanz, R. G. Fernández, and C. Vega, A potential model for the study of ices and amorphous water: TIP4P/Ice, J. Chem. Phys. 122 (2005) 234511.
"""

import numpy as np
from math import pi,sin,cos

L1 = 0.9572 / 10
L2 = 0.1577 / 10
theta=104.52 * pi/180


hy = L1*sin(theta/2)
hz = L1*cos(theta/2)
mz = L2

sites = np.array([[0.0, 0.0, mz],
                  [0.0, hy,  hz],
                  [0.0,-hy,  hz],
                  [0.0, 0.0, 0.0]])
sites -= (sites[1]+sites[2]+sites[3]*16)/18

                  
atoms = [".","H","H","O"]
labels = ["MW","HW1","HW2","OW"]
name = "ICE"


