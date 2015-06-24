import math
import numpy as np


oh = 0.09572
hangle = 104.52 * math.pi / 180 / 2
mass=18
ohz = oh * math.cos(hangle)
ohy = oh * math.sin(hangle)



sites = np.array([[0, 0,-ohz*2/mass],
                  [0, ohy,ohz*16/mass],
                  [0,-ohy,ohz*16/mass]]) # nm, OHHM

labels = ["OW","HW1","HW2"]
