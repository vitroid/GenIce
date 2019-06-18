import logging

import numpy as np
import pairlist as pl

from genice import rigid
from genice.molecules import tip4p


def hbonds(waters, cell, rotmat):
    O = np.zeros((len(waters), 3))
    H1 = np.zeros_like(O)
    H2 = np.zeros_like(O)
    for i, (com, R) in enumerate(zip(waters, rotmat)):
        a = tip4p.sites @ R + com
        O[i] = a[0]
        H1[i] = a[1]
        H2[i] = a[2]
    celli = np.linalg.inv(cell)
    O = O @ celli
    H1 = H1 @ celli
    H2 = H2 @ celli
    grid = pl.determine_grid(cell, 0.245)
    pairs = [(i, j) for i, j in pl.pairs_fine_hetero(H1, O, 0.245, cell, grid, distance=False) if i != j]
    pairs += [(i, j) for i, j in pl.pairs_fine_hetero(H2, O, 0.245, cell, grid, distance=False) if i != j]
    return pairs


class Loader():  # for analice
    def __init__(self, filename, oname="", hname="", avgspan=0):  # oname, hname and avgspan are unused.
        logger = logging.getLogger()
        logger.debug('load {0}'.format(filename))
        self.file = open(filename)
        self.bondlen = 0.3

    def load_iter(self):
        logger = logging.getLogger()
        logger.info("  Loading NX3A assuming TIP4P water.")
        while True:
            line = self.file.readline()
            if len(line) == 0:
                return
            if len(line) > 4:
                if line[:5] == "@BOX3":
                    logger.info("  @BOX3")
                    line = self.file.readline()
                    box = np.array([float(x) for x in line.split()[:3]])
                    self.cell = np.diag(box) / 10  # in nm
                    celli = np.linalg.inv(self.cell)
                    self.celltype = 'triclinic'
                elif line[:5] == "@NX3A":
                    line = self.file.readline()
                    nmol = int(line.split()[0])
                    self.waters = []
                    logger.info("  @NX3A")
                    self.rotmat = []
                    for i in range(nmol):
                        line = self.file.readline()
                        cols = line.split()[:6]
                        euler = np.array([float(x) for x in cols[3:6]])
                        self.rotmat.append(rigid.euler2rotmat(euler))
                        pos = np.array([float(x) for x in cols[:3]])
                        self.waters.append(pos / 10)  # in nm
                    self.coord = 'absolute'
                    self.pairs = hbonds(self.waters, self.cell, self.rotmat)
                    self.density = len(self.waters) / (np.linalg.det(self.cell) * 1e-21) * 18 / 6.022e23
                    yield self
