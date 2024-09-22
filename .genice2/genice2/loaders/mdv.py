import numpy as np
import re
import logging

desc = {
    "ref": {},
    "brief": "MDView file (in Angdtrom).",
    "usage": "No options available.",
}

if __name__[-4:] == "mdva":
    desc["brief"] = "MDView file (in au)."


def load_iter(file, oname, hname=None):
    logger = logging.getLogger()
    if __name__[-4:] == "mdva":
        conv = 0.052917721067  # au in nm
    else:
        conv = 0.1  # AA in nm

    natom = -100
    lineno = 0
    for line in iter(file.readline, ""):
        lineno += 1
        if lineno == 1:
            # first line
            cols = line.split()
            assert cols[0] == "#"  # yaga style
            c = [float(x) for x in cols[1:4]]
            cellmat = np.array([[c[0], 0.0, 0.0], [0.0, c[1], 0.0], [0.0, 0.0, c[2]]])
            cellmat *= conv
            celli = np.linalg.inv(cellmat)
        elif line[0] == "-":
            continue
        elif natom < 0:
            natom = int(line)
            hatoms = []
            oatoms = []
            skipped = set()
        elif natom > 0:
            natom -= 1
            cols = line.split()
            atomname = cols[0]
            # atomid = int(line[15:20])
            pos = np.array([float(x) for x in cols[1:4]]) * conv
            if re.fullmatch(oname, atomname):
                oatoms.append(pos)
            elif hname is not None and re.fullmatch(hname, atomname):
                hatoms.append(pos)
            else:
                if atomname not in skipped:
                    logger.info("Skip {0}".format(atomname))
                    skipped.add(atomname)
            # finalize the frame
            if natom == 0:
                # fractional coordinate
                oatoms = np.array(oatoms) @ celli
                if len(hatoms) > 0:
                    hatoms = np.array(hatoms) @ celli
                else:
                    hatoms = None

                yield oatoms, hatoms, cellmat
                natom = -100
                lineno = 0
