# coding: utf-8
"""
re-make python module for GenIce
"""

import numpy as np
import itertools as it
from math import floor
MAXCELL=11
torr = 1e-8

def isZero(x):
    return -torr < x < torr



def OccupyCell(nv, cells, ijk):
    if nv not in cells:
        #print("kill:",nv)
        cells.add(nv)
        for d in it.product((-2,-1,0,1,2), repeat=3):
            vv = tuple(np.dot(d,ijk)+nv)
            if d != (0,0,0):
                #print("pin", vv, d)
                #if vv in cells:
                #    print("??", vv)
                cells.add(vv)


def FindEmptyCell(cellmat, coord, ijk, labels=None):
    newcell = np.dot(ijk, cellmat)
    L1 = np.linalg.norm(newcell[0])
    L2 = np.linalg.norm(newcell[1])
    L3 = np.linalg.norm(newcell[2])
    rL = np.array([L1,L2,L3])
    # print(rL)
    ncell = 0
    cells = set()
    queue = [(0,0,0)]
    s = ""
    while len(queue) > 0:
        nv = queue.pop(0)
        if nv not in cells:
            ncell += 1
            for i, xyz in enumerate(coord):
                xxv = np.dot(xyz + nv, cellmat) # rel to abs coord
                p0v = np.dot(xxv, newcell.T)    # inner products with axes vectors
                # print(xxv,p0v,nv)
                pv  = p0v / rL / rL             # position relative to the newcell
                # print(pv,rL)
                pv -= np.floor(pv)
                # print(nv)
                label = ""
                if labels is not None:
                    label = labels[i]
                s += "{3} {0:9.4f} {1:9.4f} {2:9.4f}\n".format(pv[0],pv[1],pv[2], label)
            #print("NV:",nv)
            OccupyCell(nv, cells, ijk)
            for xi,yi,zi in it.product((-1,0,1), repeat=3):
                nei = nv[0]+xi,nv[1]+yi,nv[2]+zi
                if nei not in cells:
                    #print("nei", nei)
                    queue.append(nei)
    return ncell, s




# 2018-1-18 Still do not work with non-cubic cells.
    

ijk = np.array([[1, 1, 1], [1, -1, 0], [1, 1, -2]])
# ijk = np.array([[2, 0, 1], [0, 1, 0], [-1, 0, 2]])
# ijk = np.array([[0, 1, 0], [-1, 0, 0], [0, 0, 1]])
# ijk = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])

def hook1(lattice):
    lattice.logger.info("Hook1: Output as a python module.")
    cellmat = lattice.repcell.mat
    assert cellmat[2,0] == 0.0
    assert np.linalg.det(cellmat) > 0.0
    assert np.dot(ijk[0],ijk[1]) == 0.0
    assert np.dot(ijk[1],ijk[2]) == 0.0
    assert np.dot(ijk[2],ijk[0]) == 0.0
    lattice.logger.info("  Reshaping the unit cell.")
    lattice.logger.info("    i:{0}".format(ijk[0]))
    lattice.logger.info("    j:{0}".format(ijk[1]))
    lattice.logger.info("    k:{0}".format(ijk[2]))
    newcell = np.dot(ijk, cellmat)
    # print(newcell)
    vol = abs(np.linalg.det(ijk))
    vol = floor(vol*8192+0.5)/8192
    e1 = newcell[0].astype(float)
    e1 /= np.linalg.norm(e1)
    e3 = np.cross(newcell[0], newcell[1]).astype(float)
    e3 /= np.linalg.norm(e3)
    e2 = np.cross(e3,e1)
    R = np.array([e1,e2,e3])
    #print(R)
    RI = np.linalg.inv(R)
    #print(RI)
    #Regularization
    regcell = np.dot(newcell, RI)
    #print(regcell)
    # print(np.dot(R,R))
    #print(riv)
    #print(np.dot(riv, RI))
    #print(np.dot(rjv, RI))
    #print(np.dot(rkv, RI))
    # header
    s = ""
    s += '"""\n'
    s += "\n".join(lattice.doc) + "\n"
    s += '"""\n'
    s += "bondlen={0}\n".format(lattice.bondlen*10)
    s += "coord='relative'\n"
    if isZero(regcell[1,0]) and isZero(regcell[2,0]) and isZero(regcell[2,1]):
        s += "celltype='rect'\n"
        s += "cell='{0} {1} {2}'\n".format(regcell[0,0]*10,regcell[1,1]*10,regcell[2,2]*10)
    else:
        s += "celltype='triclinic'\n"
        s += "cell='"
        for d in range(3):
            s += "{0} {1} {2} ".format(regcell[d,0]*10,regcell[d,1]*10,regcell[d,2]*10)
        s += "'\n"
    # s += "cell='{0} {1} {2}'\n".format(ri,rj,rk)
    s += "density={0}\n".format(lattice.density)
    s += 'waters="""'+"\n"
    lattice.logger.info("  Total number of molecules: {0}".format(vol*len(lattice.reppositions)))
    
    ncell, ss = FindEmptyCell(cellmat, lattice.reppositions, ijk)
    assert vol == ncell

    #footer
    s += ss + '"""' + "\n\n"

    if lattice.cagepos is not None:
        s += 'cages="""'+"\n"
        ncell, ss = FindEmptyCell(cellmat, lattice.repcagepos, ijk, labels=lattice.repcagetype)
        s += ss + '"""'+"\n\n"
    
    print(s,end="")

    lattice.logger.info("Hook1: end.")


hooks = {1:hook1}
