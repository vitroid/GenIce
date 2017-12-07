import logging
import numpy as np

from genice.cells import rel_wrap

def Alkyl(cpos, root, cell, molname, backbone):
    """
    put a normal-alkyl group rooted at root toward cpos.
    """
    logger = logging.getLogger()
    # logger.info("  Put butyl at {0}".format(molname))
    v1abs = cell.rel2abs(rel_wrap(cpos - root))
    v1 = v1abs / np.linalg.norm(v1abs)

    origin = cell.rel2abs(root)
    CC = 0.154
    rawatoms = build(v1, v1abs*1.5/CC, backbone)

    atoms = []
    for i, atom in enumerate(rawatoms):
        atomname, pos = atom
        atompos = cell.abs_wrapf(pos*CC + origin)
        atoms.append([i, molname, atomname, atompos, 0])
    return atoms




#!/usr/bin/env python3


import numpy as np
import logging
from math import pi, sin, cos

#a sample for general alkyl group with methyls
# tree | branch: [root,branch1,branch2,branch3]
# or             [root,branch1,branch2]
# or             [root,branch1]
# or             [leaf]
# or             leaf
# 3-methylbutyl : backbone=["Ma",["Mb",["Mc",["Md","Me"]]]]]
#v1 and v2 must be given as a unit vector.
def build(direction, destination, tree, dest=None):
    """
    put a normal-alkyl group rooted at root toward cpos.
    """
    logger = logging.getLogger()
    logger.debug("  {0}".format(tree))
    # logger.info("  Put butyl at {0}".format(molname))

    if type(tree) is not list:
        tree = [tree]

    # v2 is a vector from direction to the destination
    v2 = destination - direction
    v2 /= np.linalg.norm(v2)
    v2d = np.dot(direction, v2)
    while v2d > 0.999:
        #They are inline. It is not safe to determine the orientation.
        v2 = np.random.random(3)
        v2 /= np.linalg.norm(v2)
        v2d = np.dot(direction, v2)

    v2 -= v2d * direction
    v2 /= np.linalg.norm(v2)

    #v1 is the pivot
    v1 = direction
    
    v3 = np.cross(v1, v2)     # the thild unit vector

    c = cos(120 * pi / 180)
    s = sin(120 * pi / 180)
    v4 = v2*c + v3*s   # a branch vector
    v5 = v2*c - v3*s   # another branch vector
    logger.debug("  {0} -0.5".format(np.dot(v2, v4)))
    logger.debug("  {0} -0.5".format(np.dot(v4, v5)))
    logger.debug("  {0} -0.5".format(np.dot(v5, v2)))

    c = cos(109.5 * pi / 180)
    s = sin(109.5 * pi / 180)
    v2 = -v1*c + v2*s
    v4 = -v1*c + v4*s
    v5 = -v1*c + v5*s
    logger.debug("  {0} 1/3".format(np.dot(v1, v2)))
    logger.debug("  {0} 1/3".format(np.dot(v1, v4)))
    logger.debug("  {0} 1/3".format(np.dot(v1, v5)))
    logger.debug("  {0} 1".format(np.dot(v2, v2)))
    logger.debug("  {0} 1".format(np.dot(v4, v4)))
    logger.debug("  {0} 1".format(np.dot(v5, v5)))

    #assert False
    atomname = tree[0]
    atoms = [(atomname, np.zeros(3))]
    for vec,topo in zip([v2,v4,v5], tree[1:]):
        atoms += build(vec, destination - v1, topo)
    #untranslation
    atoms = [(atom[0], atom[1] + v1) for atom in atoms]
    
    return atoms


def test():
    logging.basicConfig(level=logging.DEBUG,
                            format="%(asctime)s %(levelname)s %(message)s")
    direction = np.array([1.0, 0.0, 0.0])    #must be a unit vector
    destination = np.array([10.0, 10.0, 0.0])  #All the branches direct to this point.
    atoms = build(direction, destination, tree=["Ma",["Mb", "Mf", ["Mc", "Me", ["Md", ["A", ["B", ["C", ["D"]]]]]]]])
    # in yaplot format
    print("t 0 0 0 +")
    print("t",destination[0], destination[1], destination[2], "@")
    for atom in atoms:
        name, pos = atom
        print("t",pos[0], pos[1], pos[2], name)


if __name__ == "__main__":
    test()
