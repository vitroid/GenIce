import importlib
import logging
import pkg_resources as pr
import sys
from genice.cell import Cell
import numpy as np
from math import pi, acos

from logging import getLogger, StreamHandler, DEBUG, INFO
logger = getLogger(__name__)
handler = StreamHandler()
handler.setLevel(INFO)
logger.setLevel(INFO)
logger.addHandler(handler)
logger.propagate = False
logging.basicConfig(level=logging.INFO,
                    format="%(levelname)s %(message)s")

name = sys.argv[1]
module = importlib.import_module(name)

cell = Cell(module.cell.mat)
La = np.linalg.norm(cell.mat[0])
Lb = np.linalg.norm(cell.mat[1])
Lc = np.linalg.norm(cell.mat[2])
alpha = acos((cell.mat[1] @ cell.mat[2]) / (Lb * Lc)) * 180 / pi
beta = acos((cell.mat[2] @ cell.mat[0]) / (Lc * La)) * 180 / pi
gamma = acos((cell.mat[0] @ cell.mat[1]) / (La * Lb)) * 180 / pi

s = []
s.append("a={0}".format(La))
s.append("b={0}".format(Lb))
s.append("c={0}".format(Lc))

if alpha != 90.0:
    s.append("A={0}".format(alpha))

if beta != 90.0:
    s.append("B={0}".format(beta))

if gamma != 90.0:
    s.append("C={0}".format(gamma))

for line in open(name + ".py").readlines():
    if line.find("celltype") >= 0:
        pass
    else:
        print(line, end="")

print()
print("from genice.cell import cellvectors")
H = "cell = cellvectors("
print(H + (",\n" + " " * len(H)).join(s) + ")")
