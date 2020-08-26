desc={"ref": {"gro": "http://manual.gromacs.org/current/online/gro.html"},
      "brief": "Gromacs .gro file.",
      "usage": "No options available."
      }


import logging
import re

import numpy as np
import pairlist as pl

from genice2.cell import rel_wrap


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


def load_iter(file, oname="O", hname=None):
    logger = logging.getLogger()
    while True:
        line = readaline(file)
        if len(line) == 0:
            return
        line = readaline(file)
        if len(line) == 0:
            return
        natom = int(line)
        hatoms = []
        oatoms = []
        skipped = set()
        for i in range(natom):
            line = file.readline()
            # resid = int(line[0:5])
            # resna = line[5:10]
            atomname = line[10:15].replace(' ', '')
            # atomid = int(line[15:20])
            pos = np.array([float(x) for x in line[20:].split()[:3]])  # drop velocity
            if re.fullmatch(oname, atomname):
                oatoms.append(pos)
            elif hname is not None and re.fullmatch(hname, atomname):
                hatoms.append(pos)
            else:
                if atomname not in skipped:
                    logger.info("Skip {0}".format(atomname))
                    skipped.add(atomname)
        c = [float(x) for x in file.readline().split()]
        if len(c) == 3:
            cellmat = np.array([[c[0], 0., 0.],
                                [0., c[1], 0.],
                                [0., 0., c[2]]])
        else:
            cellmat = np.array([[c[0], c[3], c[4]],
                                [c[5], c[1], c[6]],
                                [c[7], c[8], c[2]]])
        celli = np.linalg.inv(cellmat)
        oatoms = np.array(oatoms) @ celli
        if len(hatoms) == 0:
            hatoms = None
        else:
            hatoms = np.array(hatoms) @ celli
        yield oatoms, hatoms, cellmat
