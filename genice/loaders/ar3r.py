import logging

import numpy as np


class Loader():  # for analice
    def __init__(self, filename, oname="", hname="", avgspan=0):  # oname, hname and avgspan are unused.
        logger = logging.getLogger()
        logger.debug('load {0}'.format(filename))
        self.file = open(filename)
        self.bondlen = 0.3

    def load_iter(self):
        logger = logging.getLogger()
        logger.info("  Loading AR3Rs...")
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
                elif line[:5] == "@AR3R":
                    line = self.file.readline()
                    nmol = int(line.split()[0])
                    self.waters = []
                    logger.info("  @AR3R")
                    for i in range(nmol):
                        line = self.file.readline()
                        cols = line.split()[:3]
                        pos = np.array([float(x) for x in cols[:3]])
                        self.waters.append(pos)
                    self.coord = 'relative'
                    self.rotmat = None
                    self.density = len(self.waters) / (np.linalg.det(self.cell) * 1e-21) * 18 / 6.022e23
                    yield self
