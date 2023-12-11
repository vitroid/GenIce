"""
Parallelepiped cell
"""

import numpy as np
import logging
from math import pi, sin, cos, sqrt, acos, degrees
import itertools as it


## Moved from formats/reshape.py
MAXCELL = 11
torr = 1e-8


def isZero(x):
    return -torr < x < torr


def FlagEquivCells(nv, flags, ijk):
    if nv not in flags:
        for d in it.product((-2, -1, 0, 1, 2), repeat=3):
            vv = tuple(d @ ijk + nv)
            flags.add(vv)


def FindEmptyCells(cellmat, ijk, relpositions, labels=None):
    newcell = ijk @ cellmat
    newcelli = np.linalg.inv(newcell)
    L1 = np.linalg.norm(newcell[0])
    L2 = np.linalg.norm(newcell[1])
    L3 = np.linalg.norm(newcell[2])
    rL = np.array([L1, L2, L3])
    # print(rL)
    ncell = 0
    flags = set()
    queue = [(0, 0, 0)]
    s = ""
    while len(queue) > 0:
        nv = queue.pop(0)
        if nv not in flags:
            ncell += 1
            # Place atoms
            for i, xyz in enumerate(relpositions):
                # rel to abs
                xxv = (xyz + nv) @ cellmat
                # inner products with axis vector of the newcell
                pv = xxv @ newcelli
                # print(pv,rL)
                pv -= np.floor(pv)
                # print(nv)
                label = ""
                if labels is not None:
                    label = labels[i]
                s += "{3} {0:9.4f} {1:9.4f} {2:9.4f}\n".format(
                    pv[0], pv[1], pv[2], label
                )
            # print("NV:",nv)
            FlagEquivCells(nv, flags, ijk)
            for xi, yi, zi in it.product((-1, 0, 1), repeat=3):
                nei = nv[0] + xi, nv[1] + yi, nv[2] + zi
                if nei not in flags:
                    # print("nei", nei)
                    queue.append(nei)
    return ncell, s


def rel_wrap(relvec):
    return relvec - np.floor(relvec + 0.5)


def rel_wrapf(relvec):
    return relvec - np.floor(relvec)


def sincos(A):
    cA = cos(A)
    if abs(cA) < 1e-16:
        return 1.0, 0.0
    return sin(A), cA


# moved from GenIce/cif/CIF.py
def cellvectors(a, b, c, A=90, B=90, C=90):
    """
    Generate cell vectors from a,b,c and alpha, beta, gamma.
    """
    # probably same as six2nine in rigid.py
    logger = logging.getLogger()
    A *= pi / 180
    B *= pi / 180
    C *= pi / 180
    sA, cA = sincos(A)
    sB, cB = sincos(B)
    sC, cC = sincos(C)
    ea = np.array([1.0, 0.0, 0.0])
    eb = np.array([cC, sC, 0])
    # ec.ea = ecx = cos(B)
    # ec.eb = ecx*ebx + ecy*eby = cos(A)
    ecx = cB
    ecy = (cA - ecx * eb[0]) / eb[1]
    ecz = sqrt(1 - ecx**2 - ecy**2)
    ec = np.array([ecx, ecy, ecz])
    return np.vstack([ea * a, eb * b, ec * c])


def cellshape(cellmat):
    """
    From cell matrix to a, b, c, alpha, beta, and gamma.
    """
    a = np.linalg.norm(cellmat[0])
    b = np.linalg.norm(cellmat[1])
    c = np.linalg.norm(cellmat[2])
    alpha = degrees(acos((cellmat[1] @ cellmat[2]) / (b * c)))
    beta = degrees(acos((cellmat[2] @ cellmat[0]) / (c * a)))
    gamma = degrees(acos((cellmat[0] @ cellmat[1]) / (a * b)))
    return a, b, c, alpha, beta, gamma


class Cell:
    def __init__(self, desc=None):
        self.mat = np.zeros([3, 3])
        self.inv = None
        self.parse(desc)

    def abs2rel(self, absvecs):
        return absvecs @ self.inv

    def rel2abs(self, relvec):
        return relvec @ self.mat

    def abs_wrap(self, absvec):
        return self.rel2abs(rel_wrap(self.abs2rel(absvec)))

    def abs_wrapf(self, absvec):
        return self.rel2abs(rel_wrapf(self.abs2rel(absvec)))

    def volume(self):
        return np.linalg.det(self.mat)

    def scale(self, x):
        self.mat *= x
        self.inv = np.linalg.inv(self.mat)

    def scale2(self, x):
        for d in range(3):
            self.mat[d, :] = self.mat[d, :] * x[d]
        self.inv = np.linalg.inv(self.mat)

    def parse(self, mat):
        logger = logging.getLogger()
        logger.debug(("MAT:", mat))
        self.mat = mat.copy()
        La = np.linalg.norm(self.mat[0])
        Lb = np.linalg.norm(self.mat[1])
        Lc = np.linalg.norm(self.mat[2])
        alpha = acos((self.mat[1] @ self.mat[2]) / (Lb * Lc)) * 180 / pi
        beta = acos((self.mat[2] @ self.mat[0]) / (Lc * La)) * 180 / pi
        gamma = acos((self.mat[0] @ self.mat[1]) / (La * Lb)) * 180 / pi
        logging.info("Cell dimension:")
        logging.info("  a = {0}".format(La))
        logging.info("  b = {0}".format(Lb))
        logging.info("  c = {0}".format(Lc))
        logging.info("  A = {0}".format(alpha))
        logging.info("  B = {0}".format(beta))
        logging.info("  C = {0}".format(gamma))
        self.inv = np.linalg.inv(self.mat)

    def serialize_BOX9(self):
        s = "@BOX9\n"
        for d in range(3):
            s += "{0:.4f} {1:.4f} {2:.4f}\n".format(*self.mat[d] * 10)  # AA
        return s

    def shape(self):
        return cellshape(self.mat)
