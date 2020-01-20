from logging import getLogger

import numpy as np

from genice import rigid
from genice.molecules import tip4p


def load_iter(file, **kwargs):
    logger = getLogger()
    logger.info("  Loading NX4A assuming TIP4P water.")
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
            elif line[:5] == "@NX4A":
                line = file.readline()
                nmol = int(line.split()[0])
                oatoms = []
                hatoms = []
                logger.info("  @NX4A")
                for i in range(nmol):
                    line = file.readline()
                    cols = line.split()[:7]
                    pos = np.array([float(x) for x in cols[:3]]) / 10
                    # oatoms.append(pos / 10)  # in nm
                    quat = np.array([float(x) for x in cols[3:7]])
                    rotmat = rigid.quat2rotmat(quat)
                    a = tip4p.sites @ rotmat + pos
                    oatoms.append(a[0])
                    hatoms.append(a[1])
                    hatoms.append(a[2])
                oatoms = np.array(oatoms) @ celli
                hatoms = np.array(hatoms) @ celli
                yield oatoms, hatoms, cellmat
