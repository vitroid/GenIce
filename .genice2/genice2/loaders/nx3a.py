import logging

import numpy as np

from genice2 import rigid
import genice2.molecules.tip4p


def load_iter(file, **kwargs):
    logger = logging.getLogger()
    logger.info("  Loading NX3A assuming TIP4P water.")
    tip4p = genice2.molecules.tip4p.Molecule()
    name, labels, sites = tip4p.get()
    while True:
        line = file.readline()
        if len(line) == 0:
            return
        if len(line) > 4:
            if line[:5] == "@BOX3":
                logger.info("  @BOX3")
                line = file.readline()
                box = np.array([float(x) for x in line.split()[:3]])
                cellmat = np.diag(box) / 10  # in nm
                celli = np.linalg.inv(cellmat)
            elif line[:5] == "@NX3A":
                line = file.readline()
                nmol = int(line.split()[0])
                oatoms = []
                hatoms = []
                logger.info("  @NX3A")
                for i in range(nmol):
                    line = file.readline()
                    cols = line.split()[:6]
                    pos = np.array([float(x) for x in cols[:3]]) / 10
                    # oatoms.append(pos / 10)  # in nm
                    euler = np.array([float(x) for x in cols[3:6]])
                    rotmat = rigid.euler2rotmat(euler)
                    a = sites @ rotmat + pos
                    oatoms.append(a[0])
                    hatoms.append(a[1])
                    hatoms.append(a[2])
                oatoms = np.array(oatoms) @ celli
                hatoms = np.array(hatoms) @ celli
                yield oatoms, hatoms, cellmat
