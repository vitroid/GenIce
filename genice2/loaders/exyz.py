from genice2.cell import rel_wrap
import numpy as np
import re
from logging import getLogger
desc = {
    "ref": {
        "exyz": "https://open-babel.readthedocs.io/en/latest/FileFormats/Extended_XYZ_cartesian_coordinates_format.html"},
    "brief": "Extended XYZ format.",
    "usage": "No options available."}


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
    logger = getLogger()
    while True:
        # first line: number of atoms
        line = readaline(file)
        if len(line) == 0:  # EOF
            break
        # second line: %PBC
        line = readaline(file)
        if len(line) == 0:  # EOF
            break
        assert line[:4] == "%PBC"
        hatoms = []
        oatoms = []
        skipped = set()
        while True:
            line = file.readline()
            cols = line.split()
            if len(cols) == 0:
                break
            # here we treat the data as free format, although
            # each column width is fixed according to the refernce.
            atomname = cols[0]
            # atomid = int(line[15:20])
            pos = np.array([float(x) / 10 for x in cols[1:4]])  # in nm
            if re.fullmatch(oname, atomname):
                oatoms.append(pos)
            elif hname is not None and re.fullmatch(hname, atomname):
                hatoms.append(pos)
            else:
                if atomname not in skipped:
                    logger.info("Skip {0}".format(atomname))
                    skipped.add(atomname)
        line = readaline(file)
        assert line[:7] == "Vector1"
        v1 = [float(x) for x in line.split()[1:4]]
        line = readaline(file)
        assert line[:7] == "Vector2"
        v2 = [float(x) for x in line.split()[1:4]]
        line = readaline(file)
        assert line[:7] == "Vector3"
        v3 = [float(x) for x in line.split()[1:4]]
        line = readaline(file)
        assert line[:6] == "Offset"
        offs = [float(x) for x in line.split()[1:4]]
        # The meaning of "Offset" is anbiguous.
        # We will ignore it for now.

        cellmat = np.array([v1, v2, v3]) / 10  # in nm
        celli = np.linalg.inv(cellmat)
        # fractional coordinate
        oatoms = np.array(oatoms) @ celli
        if len(hatoms) == 0:
            hatoms = None
        else:
            # fractional coordinate
            hatoms = np.array(hatoms) @ celli
        yield oatoms, hatoms, cellmat
    return


"""
4
%PBC
   C        0.00000        1.40272        0.00000
   H        0.00000        2.49029        0.00000
   C       -1.21479        0.70136        0.00000
   H       -2.15666        1.24515        0.00000

Vector1    2.445200    0.000000    0.000000
Vector2    0.000000    1.000000    0.000000
Vector3    0.000000    0.000000    1.000000
Offset     0.000000    0.000000    0.000000
"""
