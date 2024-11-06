from genice2.cell import rel_wrap
import numpy as np
import re
import logging

desc = {
    "ref": {"gro": "http://manual.gromacs.org/current/online/gro.html"},
    "brief": "Gromacs .gro file.",
    "usage": "No options available.",
}


def readaline(file):
    """
    read a non-comment line.
    """
    for line in iter(file.readline(), ""):
        if line[0] != "#":
            return line


def load_iter(file, oname="O", hname=None):
    logger = logging.getLogger()

    natom = -100
    lineno = 0
    for line in iter(file.readline, ""):
        # skip whenever the line starts from "#"
        if len(line) > 0 and line[0] == "#":
            continue
        lineno += 1
        if lineno == 1:
            logging.debug(f"COMMENT {line.rstrip()}")
            # the first line is a comment
            continue
        elif lineno == 2:
            # the second line is th enumber of atoms
            natom = int(line)
            logger.debug(f"NATOM {natom}")
            hatoms = []
            oatoms = []
            skipped = set()
        elif natom > 0:
            # count down the natom to read atoms
            natom -= 1
            logger.debug(f"NATOM {natom}")
            atomname = line[10:15].replace(" ", "")
            pos = np.array([float(x) for x in line[20:].split()[:3]])  # drop velocity
            if re.fullmatch(oname, atomname):
                oatoms.append(pos)
            elif hname is not None and re.fullmatch(hname, atomname):
                hatoms.append(pos)
            else:
                skipped.add(atomname)
                if atomname not in skipped:
                    logger.info("Skip {0}".format(atomname))

        elif natom == 0:
            # the last line is the cell dimension
            c = [float(x) for x in line.split()]
            if len(c) == 3:
                cellmat = np.diag(c)
            else:
                cellmat = np.array(
                    [[c[0], c[3], c[4]], [c[5], c[1], c[6]], [c[7], c[8], c[2]]]
                )
            celli = np.linalg.inv(cellmat)
            oatoms = np.array(oatoms) @ celli
            if len(hatoms) == 0:
                hatoms = None
            else:
                hatoms = np.array(hatoms) @ celli
            yield oatoms, hatoms, cellmat

            # reset everything
            natom = -100
            lineno = 0
