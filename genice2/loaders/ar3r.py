import logging

import numpy as np


def load_iter(file, **kwargs):
    logger = logging.getLogger()
    logger.info("  Loading AR3Rs...")
    cellmat = None
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
            elif line[:5] == "@AR3R":
                line = file.readline()
                nmol = int(line.split()[0])
                oatoms = []
                logger.info("  @AR3R")
                for i in range(nmol):
                    line = file.readline()
                    cols = line.split()[:3]
                    pos = np.array([float(x) for x in cols[:3]])
                    oatoms.append(pos)
                hatoms = None
                oatoms = np.array(oatoms)
                yield oatoms, hatoms, cellmat
