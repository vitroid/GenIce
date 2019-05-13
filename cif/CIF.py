#!/usr/bin/python

import numpy as np
from math import *
import pairlist as pl
from collections import defaultdict
import logging

from logging import getLogger, StreamHandler, DEBUG, INFO
logger = getLogger(__name__)
handler = StreamHandler()
handler.setLevel(INFO)
logger.setLevel(INFO)
logger.addHandler(handler)
logger.propagate = False



def fullatoms(atomd, sops):
    global x,y,z
    # the variables to be evaluated by eval() must be global (?)
    full = []
    names = []
    for name, pos in atomd.items():
        x,y,z = pos
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
                    tooclose = True
                    break
            if not tooclose:
                yield name, p
                full.append(p)
                names.append(name)


def cellvectors(a,b,c,A=90,B=90,C=90):
    A *= pi/180
    B *= pi/180
    C *= pi/180
    ea = np.array([1.0, 0.0, 0.0])
    eb = np.array([cos(C), sin(C), 0])
    # ec.ea = ecx = cos(B)
    # ec.eb = ecx*ebx + ecy*eby = cos(A)
    ecx = cos(B)
    ecy = (cos(A) - ecx*eb[0]) / eb[1]
    ecz = sqrt(1-ecx**2-ecy**2)
    ec = np.array([ecx, ecy, ecz])
    logger.debug((cos(A), np.dot(eb, ec)))
    logger.debug((cos(B), np.dot(ec, ea)))
    logger.debug((cos(C), np.dot(ea, eb)))
    return np.vstack([ea*a, eb*b, ec*c])


def atomdic(atoms):
    atomd = dict()
    for atom in atoms.split("\n"):
        cols = atom.split()
        if len(cols) == 0:
            continue
        # remove error digits from 1-3 cols.
        xyz = np.array([float(col[:col.find("(")]) for col in cols[1:]])
        name = cols[0]
        atomd[name] = xyz
    return atomd


def symmetry_operators(symops):
    sops = []
    for symop in symops.split("\n"):
        cols = symop.split()
        if len(cols) == 0:
            continue
        sops.append(cols)
    return sops


def gromacs(cell, atomd, sops):
    oxygens = []
    hydrogens = []
    for name, pos in fullatoms(atomd, sops):
        if name[0] in ("O",):
            oxygens.append(pos)
        elif name[0] in "DH":
            hydrogens.append(pos)

    oxygens = np.array(oxygens)
    hydrogens = np.array(hydrogens)
    # print(oxygens.shape, hydrogens.shape)

    oh = defaultdict(list)
    grid = pl.determine_grid(cell, 0.12)
    for i,j in pl.pairs_fine_hetero(oxygens, hydrogens, 0.12, cell, grid, distance=False):
        oh[i].append(j)


    class Lattice():
        def __init__(self):
            pass

    import genice.cells as ce

    lattice = Lattice()
    lattice.logger = logging.getLogger()
    lattice.repcell = ce.Cell(cell, celltype='triclinic')
    lattice.doc = ""

    lattice.atoms = []
    for i,O in enumerate(oh):
        H0, H1 = oh[O]
        lattice.atoms.append([i+1, "SOL", "O",  oxygens[O] @ cell,    0])
        lattice.atoms.append([i+1, "SOL", "H1", hydrogens[H0] @ cell, 0])
        lattice.atoms.append([i+1, "SOL", "H2", hydrogens[H1] @ cell, 0])

    from genice.formats.gromacs import hook7

    hook7(lattice)


