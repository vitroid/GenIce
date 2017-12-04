#!/usr/bin/env python3


import numpy as np
import logging
from math import pi, sin, cos

#a sample for general alkyl group with methyls
# 3-methylbutyl : backbone=[["Ma",[],[]],["Mb",[],[]],["Mc",[["Me",[],[]]],[]],["Md",[],[]]]
#v1 and v2 must be given as a unit vector.
def alkyl(v1, v2, backbone):
    """
    put a normal-alkyl group rooted at root toward cpos.
    """
    logger = logging.getLogger()
    # logger.info("  Put butyl at {0}".format(molname))
    v3 = np.cross(v1, v2)     # the thild unit vector

    c = cos(109.5 / 2 * pi / 180)
    s = (1.0 - c**2)**0.5
    v4 = v2*c + v3*s   # a branch vector
    v5 = v2*c - v3*s   # another branch vector

    logger.debug("  {0} {1} {2}".format(
        np.dot(v1, v1), np.dot(v2, v2), np.dot(v1, v2)))

    atoms = []
    for i, meth in enumerate(backbone):
        center, b1, b2 = meth
        x = s * (i + 1)
        y = ((i + 1) % 2) * c
        cpos = x * v1 + y * v2
        atoms.append([center, cpos])

        if len(b1) > 0:
            newv1 = (v4-v5)/2
            newv2 = v4 - newv1
            newv1 /= np.linalg.norm(newv1)
            newv2 /= np.linalg.norm(newv2)
            adatoms = alkyl(newv1, newv2, b1)
            for atom in adatoms:
                atoms.append((atom[0],atom[1]+cpos))
        if len(b2) > 0:
            newv1 = (v5-v4)/2
            newv2 = v5 - newv1
            newv1 /= np.linalg.norm(newv1)
            newv2 /= np.linalg.norm(newv2)
            adatoms = alkyl(newv1, newv2, b2)
            for atom in adatoms:
                atoms.append((atom[0],atom[1]+cpos))
        #Alternated
        v4 = -v4
        v5 = -v5
    return atoms


def test():
    v1 = np.array([1.0, 0.0, 0.0])
    v2 = np.array([0.0, 1.0, 0.0])
    atoms = alkyl(v1, v2, backbone=[["Ma",[],[]],["Mb",[["Mf",[],[]]],[]],["Mc",[["Me",[],[]]],[]],["Md",[],[]]])
    # in yaplot format
    for atom in atoms:
        name, pos = atom
        print("t",pos[0], pos[1], pos[2], name)


if __name__ == "__main__":
    test()
