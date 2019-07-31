import logging
import re

import numpy as np
import pairlist as pl

from genice.cell import rel_wrap, Cell


def readaline(file):
    """
    read a non-comment line.
    """
    while True:
        line = file.readline()
        if len(line) == 0:
            return line
        if line[0] != "#":
            return line


class Loader():  # for analice
    def __init__(self, filename, oname, hname):
        self.file = open(filename)
        self.oname = oname
        self.hname = hname

    def load_iter(self):
        logger = logging.getLogger()
        while True:
            line = readaline(self.file)
            if len(line) == 0:
                return
            line = readaline(self.file)
            if len(line) == 0:
                return
            natom = int(line)
            hatoms = []
            oatoms = []
            self.waters = []
            skipped = set()
            for i in range(natom):
                line = self.file.readline()
                # resid = int(line[0:5])
                # resna = line[5:10]
                atomname = line[10:15].replace(' ', '')
                # atomid = int(line[15:20])
                pos = np.array([float(x) for x in line[20:].split()[:3]])  # drop velocity
                if atomname == self.oname:
                    oatoms.append(pos)
                elif self.hname is not None and re.fullmatch(self.hname, atomname):
                    hatoms.append(pos)
                else:
                    if atomname not in skipped:
                        logger.info("Skip {0}".format(atomname))
                        skipped.add(atomname)
            c = [float(x) for x in self.file.readline().split()]
            if len(c) == 3:
                self.cell = np.array([[c[0], 0., 0.],
                                      [0., c[1], 0.],
                                      [0., 0., c[2]]])
            else:
                self.cell = np.array([[c[0], c[3], c[4]],
                                      [c[5], c[1], c[6]],
                                      [c[7], c[8], c[2]]])
            self.celltype = 'triclinic'
            self.coord = 'absolute'
            self.density = len(oatoms) / (np.linalg.det(self.cell) * 1e-21) * 18 / 6.022e23
            celli = np.linalg.inv(self.cell)
            ro = np.array([np.dot(x, celli) for x in oatoms])
            if len(hatoms) > 0:
                rh = np.array([np.dot(x, celli) for x in hatoms])
            self.waters = np.dot(ro, self.cell)  # abs pos
            if len(hatoms) > 0:
                self.rotmat = []
                for i in range(len(self.waters)):
                    rdh0 = rel_wrap(rh[i * 2] - ro[i])
                    rdh1 = rel_wrap(rh[i * 2 + 1] - ro[i])
                    o = np.dot(ro, self.cell)
                    dh0 = np.dot(rdh0, self.cell)
                    dh1 = np.dot(rdh1, self.cell)
                    y = dh0 - dh1
                    y /= np.linalg.norm(y)
                    z = dh0 + dh1
                    z /= np.linalg.norm(z)
                    x = np.cross(y, z)
                    self.rotmat.append(np.vstack([x, y, z]))
                    # 重心位置を補正。
                    self.waters[i] += (dh0 + dh1) * 1. / 18.
                grid = pl.determine_grid(self.cell, 0.245)
                # remove intramolecular OHs
                # 水素結合は原子の平均位置で定義している。
                self.pairs = []
                logger.debug("  Make pair list.")
                for o, h in pl.pairs_fine_hetero(ro, rh, 0.245, self.cell, grid, distance=False):
                    if not (h == o * 2 or h == o * 2 + 1):
                        # hとoは別の分子の上にあって近い。
                        # register a new intermolecular pair
                        self.pairs.append((h // 2, o))
                logger.debug("  # of pairs: {0} {1}".format(len(self.pairs), len(self.waters)))
            yield self
