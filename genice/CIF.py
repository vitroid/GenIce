#!/usr/bin/python
"""
Helpers for lattice plugin.

Generate a lattice from CIF-like data.
"""


import numpy as np
from math import *
import pairlist as pl
from collections import defaultdict

from logging import getLogger, StreamHandler, DEBUG, INFO
logger = getLogger(__name__)
handler = StreamHandler()
handler.setLevel(DEBUG)
logger.setLevel(DEBUG)
logger.addHandler(handler)
logger.propagate = False



def fullatoms(atomd, sops):
    global x,y,z
    # the variables to be evaluated by eval() must be global (?)
    full = []
    names = []
    for name, pos in atomd.items():
        x,y,z = pos
        X,Y,Z = x,y,z # alias
        for sop in sops:
            # print(x,sop)
            p = np.array([eval(s) for s in sop])
            p -= np.floor(p)
            tooclose = False
            for f in full:
                d = f - p
                d -= np.floor(d+0.5)
                L2 = np.dot(d,d)
                if L2 < 0.0001:
                    # print(f,p)
                    logger.debug("Too close: {0} {1}".format(f,p))
                    tooclose = True
                    break
            if not tooclose:
                yield name, p
                full.append(p)
                names.append(name)




def atomdic(atoms):
    atomd = dict()
    for atom in atoms.split("\n"):
        cols = atom.split()
        if len(cols) == 0:
            continue
        # remove error digits from 1-3 cols.
        xyz = []
        for col in cols[1:4]:
            p = col.find("(")
            if p >= 0:
                col = col[:p]
            xyz.append(float(col))
        name = cols[0]
        atomd[name] = xyz
    return atomd


def symmetry_operators(symops):
    sops = []
    for symop in symops.split("\n"):
        cols = symop.split()
        if len(cols) == 0:
            continue
        ops = []
        for col in cols:
            if col[-1] == ",":
                col = col[:-1]
            ops.append(col)
        sops.append(ops)
    return sops


def waters_and_pairs(cell, atomd, sops, rep=(1,1,1)):
    oxygens = []
    hydrogens = []
    for name, pos in fullatoms(atomd, sops):
        if name[0] in ("O",):
            oxygens.append(pos)
        elif name[0] in "DH":
            hydrogens.append(pos)

    assert len(oxygens)*2 == len(hydrogens) or len(hydrogens) == 0, "{0}:{1}".format(len(oxygens)*2,len(hydrogens))

    cell *= np.array(rep)
    
    oo = [[o[0]+x, o[1]+y, o[2]+z]
           for o in oxygens
           for x in range(rep[0])
           for y in range(rep[1])
           for z in range(rep[2])]
    oxygens = np.array(oo)
    oxygens /= np.array(rep)

    if len(hydrogens) == 0:
        return oxygens, None

    hh = [[h[0]+x, h[1]+y, h[2]+z]
           for h in hydrogens
           for x in range(rep[0])
           for y in range(rep[1])
           for z in range(rep[2])]
    hydrogens = np.array(hh)
    hydrogens /= np.array(rep)

    logger.debug([oxygens.shape, hydrogens.shape])

    oh = defaultdict(list)
    parent = dict()
    grid = pl.determine_grid(cell, 0.15)
    # find covalent OH bonds
    for i,j in pl.pairs_fine_hetero(oxygens, hydrogens, 0.15, cell, grid, distance=False):
        oh[i].append(j)
        parent[j] = i
    for i in oh:
        logger.debug((i,oh[i]))
    
    grid = pl.determine_grid(cell, 0.20)
    pairs = []
    # find HBs
    for i,j in pl.pairs_fine_hetero(oxygens, hydrogens, 0.20, cell, grid, distance=False):
        if j not in oh[i]:
            # H is on a different water molecule
            p = parent[j]
            pairs.append([p, i])

    logger.debug(pairs)

    waters = []
    for i in range(len(oh)):
        logger.debug((i, oh[i]))
        j,k = oh[i]
        com = (oxygens[i]*16 + hydrogens[j] + hydrogens[k])/18
        waters.append(com)

    return waters, pairs
