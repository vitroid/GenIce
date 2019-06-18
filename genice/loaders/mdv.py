import logging

import numpy as np


class Loader():  # for analice
    """
    Loader for mdview file.
    Length unit is fixed to atomic unit.
    """

    def __init__(self, filename, oname="O", hname="H", avgspan=0):  # avgspan is not used
        self.file = open(filename)
        self.oname = oname
        self.hname = hname

    def load_iter(self):
        logger = logging.getLogger()
        au = 0.052917721067  # nm
        while True:
            line = self.file.readline()  # 1st line:comment
            if len(line) == 0:
                return
            cols = line.split()
            assert cols[0] == '#'  # yaga style
            c = [float(x) for x in cols[1:4]]
            self.cell = np.array([[c[0], 0., 0.],
                                  [0., c[1], 0.],
                                  [0., 0., c[2]]])
            self.cell *= au
            line = self.file.readline()
            natom = int(line)
            hatoms = []
            self.waters = []
            skipped = set()
            for i in range(natom):
                line = self.file.readline()
                cols = line.split()
                atomname = cols[0]
                # atomid = int(line[15:20])
                pos = np.array([float(x) for x in cols[1:4]]) * au
                if atomname == self.oname:
                    self.waters.append(pos)
                elif self.hname is not None and re.fullmatch(self.hname, atomname):
                    hatoms.append(pos)
                else:
                    if atomname not in skipped:
                        logger.info("Skip {0}".format(atomname))
                        skipped.add(atomname)
            self.celltype = 'triclinic'
            self.coord = 'absolute'
            self.density = len(self.waters) / (np.linalg.det(self.cell) * 1e-21) * 18 / 6.022e23

            if len(hatoms) > 0:
                celli = np.linalg.inv(self.cell)
                # relative coord
                rh = np.array([np.dot(x, celli) for x in hatoms])
                ro = np.array([np.dot(x, celli) for x in self.waters])
                # rotmatrices for analice
                self.rotmat = []
                for i in range(len(self.waters)):
                    o = self.waters[i]
                    h0, h1 = hatoms[i * 2:i * 2 + 2]
                    h0 -= o
                    h1 -= o
                    y = h1 - h0
                    y /= np.linalg.norm(y)
                    z = h0 + h1
                    z /= np.linalg.norm(z)
                    x = np.cross(y, z)
                    self.rotmat.append(np.vstack([x, y, z]))
                grid = pl.determine_grid(self.cell, 0.245)
                # remove intramolecular OHs
                self.pairs = []
                logger.debug("  Make pair list.")
                for o, h in pl.pairs_fine_hetero(ro, rh, 0.245, self.cell, grid, distance=False):
                    if h == o * 2 or h == o * 2 + 1:
                        # adjust oxygen positions
                        dh = rh[h] - ro[o]
                        dh -= np.floor(dh + 0.5)
                        self.waters[o] += np.dot(dh, self.cell) * 1. / 16.
                    else:
                        # register a new intermolecular pair
                        self.pairs.append((h // 2, o))
                logger.debug("  # of pairs: {0} {1}".format(len(self.pairs), len(self.waters)))
            yield self
