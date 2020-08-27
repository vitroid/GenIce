#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Helpers for lattice plugin.

Generate a lattice from CIF-like data.
"""


import numpy as np
from math import *
import pairlist as pl
from collections import defaultdict

from logging import getLogger



def fullatoms(atomd, sops):
    global x,y,z
    logger = getLogger()
    # the variables to be evaluated by eval() must be global (?)
    atoms = []
    for sop in sops:
        for name, pos in atomd.items():
            x,y,z = pos
            # print(x,sop)
            p = sop(x,y,z)
            tooclose = False
            for n,f in atoms:
                d = f - p
                d -= np.floor(d+0.5)
                L2 = np.dot(d,d)
                if L2 < 0.0001:
                    # print(f,p)
                    logger.debug("Too close: {0} {1}".format(f,p))
                    tooclose = True
                    break
            if not tooclose:
                # print(p, (x,y,z), sop)
                atoms.append((name, p))
    return atoms




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
    """
    Array of symmetry operations.
    """
    def _wrap(x):
        return x - np.floor(x+0.5)

    symops = symops.translate(str.maketrans("XYZ", "xyz"))
    #symfuncs = []
    for symop in symops.split("\n"):
        cols = symop.split()
        if len(cols) == 0:
            continue
        ops = []
        for col in cols:
            if col[-1] == ",":
                col = col[:-1]
            ops.append(col)
        #symfuncs.append(lambda x,y,z: _wrap(np.array([eval(op) for op in ops])))
        yield lambda x,y,z: _wrap(np.array([eval(op) for op in ops]))
    #return symfuncs


def waters_and_pairs(cell, atomd, sops, rep=(1,1,1)):

    logger = getLogger()

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

    oh = defaultdict(list)
    parent = dict()
    # find covalent OH bonds
    for i,j in pl.pairs_iter(oxygens,
                             rc=0.15,
                             cell=cell,
                             pos2=hydrogens,
                             distance=False):
        oh[i].append(j)
        parent[j] = i
    pairs = []
    # find HBs
    for i,j in pl.pairs_iter(oxygens,
                             rc=0.20, cell=cell,
                             pos2=hydrogens,
                             distance=False):
        if j not in oh[i]:
            # H is on a different water molecule
            p = parent[j]
            pairs.append([p, i])

    logger.debug(pairs)

    waters = []
    for i in range(len(oh)):
        logger.debug((i, oh[i]))
        j,k = oh[i]
        dhj = hydrogens[j] - oxygens[i]
        dhk = hydrogens[k] - oxygens[i]
        dhj -= np.floor(dhj+0.5)
        dhk -= np.floor(dhk+0.5)
        com = oxygens[i] + (dhj + dhk)/18
        waters.append(com)

    return waters, pairs
