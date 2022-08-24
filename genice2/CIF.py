#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Helpers for lattice plugin.

Generate a lattice from CIF-like data.
"""


from collections import defaultdict
from logging import getLogger

import numpy as np
import pairlist as pl

# from math import



def fullatoms(atomd, sops):
    global x, y, z
    logger = getLogger()
    # the variables to be evaluated by eval() must be global (?)
    atoms = []
    for sop in sops:
        for name, pos in atomd.items():
            x, y, z = pos
            # print(x,sop)
            p = sop(x, y, z)
            tooclose = False
            for n, f in atoms:
                d = f - p
                d -= np.floor(d + 0.5)
                L2 = d @ d
                if L2 < 0.0001:
                    # print(f,p)
                    logger.debug("Too close: {0} {1}".format(f, p))
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


from typing import Iterable


def symmetry_operators(symops:str, offsets:Iterable=[("+0", "+0", "+0")]):
    """Generator of the symmetry operations

    Args:
        symops (str): The table of symmetry operators, e.g. http://img.chem.ucl.ac.uk/sgp/large/036az3.htm
        offset (Iterable, optional): Offset to the origin. Defaults to ["+0", "+0", "+0"].
    """
    def _wrap(x):
        return x
        return x - np.floor(x + 0.5)

    def eval_(v, x,y,z):
        # print(v)
        return eval(v)

    symops = symops.translate(str.maketrans("XYZ", "xyz"))
    #symfuncs = []
    for symop in symops.split("\n"):
        cols = symop.split(",")
        if len(cols) <= 1:
            cols = symop.split()
            if len(cols) <= 1:
                continue
        for offset in offsets:
            ops = []
            for col, offs in zip(cols, offset):
                ops.append(col+offs)
            yield lambda x, y, z: _wrap(np.array([eval_(op, x, y, z) for op in ops]))
    # return symfuncs


def waters_and_pairs(cell, atomd, sops, rep=(1, 1, 1)):

    logger = getLogger()

    oxygens = []
    hydrogens = []
    for name, pos in fullatoms(atomd, sops):
        if name[0] in ("O",):
            oxygens.append(pos)
        elif name[0] in "DH":
            hydrogens.append(pos)

    assert len(oxygens) * 2 == len(hydrogens) or len(
        hydrogens) == 0, "{0}:{1}".format(len(oxygens) * 2, len(hydrogens))

    cell *= np.array(rep)

    oo = [[o[0] + x, o[1] + y, o[2] + z]
        for o in oxygens
        for x in range(rep[0])
        for y in range(rep[1])
        for z in range(rep[2])]
    oxygens = np.array(oo)
    oxygens /= np.array(rep)

    if len(hydrogens) == 0:
        return oxygens, None

    hh = [[h[0] + x, h[1] + y, h[2] + z]
        for h in hydrogens
        for x in range(rep[0])
        for y in range(rep[1])
        for z in range(rep[2])]
    hydrogens = np.array(hh)
    hydrogens /= np.array(rep)

    oh = defaultdict(list)
    parent = dict()
    # find covalent OH bonds
    for i, j in pl.pairs_iter(oxygens,
                            maxdist=0.15, # nm
                            cell=cell,
                            pos2=hydrogens,
                            distance=False):
        oh[i].append(j)
        parent[j] = i
    pairs = []
    # find HBs
    for i, j in pl.pairs_iter(oxygens,
                            maxdist=0.20, # nm
                            cell=cell,
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
        j, k = oh[i]
        dhj = hydrogens[j] - oxygens[i]
        dhk = hydrogens[k] - oxygens[i]
        dhj -= np.floor(dhj + 0.5)
        dhk -= np.floor(dhk + 0.5)
        com = oxygens[i] + (dhj + dhk) / 18
        waters.append(com)

    return waters, pairs


# Candidates for ice XI
# {'P12_11', 'P2_12_12_1', 'Pca2_1', 'P1c1', 'Cmc2_1', 'C1c1', 'P12', 'Pbn2_1', 'Pna2_1'}
def operations(spaceg, origin=0):
    symops = None
    if spaceg == "P12_11" or spaceg == "4":
        symops = symmetry_operators("""
            x,            y,            z
            -x,          1/2+y,         -z
        """)
    elif spaceg == "P2_12_12_1" or spaceg == "19":
        # http://img.chem.ucl.ac.uk/sgp/large/019a.htm
        if origin == 0:
            symops = symmetry_operators("""
                x,            y,            z
                1/2+x,        1/2-y,         -z
                -x,          1/2+y,        1/2-z
                1/2-x,         -y,          1/2+z
            """)
    elif spaceg == "Pca2_1" or spaceg == 29:
        symops = symmetry_operators("""
            x,            y,            z
            1/2-x,          y,          1/2+z
            1/2+x,         -y,            z
            -x,           -y,          1/2+z
        """)
    elif spaceg == "P1c1" or spaceg == "7":
        symops =symmetry_operators("""
            x,            y,            z
            x,           -y,          1/2+z
        """)
    elif spaceg == "Cmc2_1" or spaceg == "36":
        symops = symmetry_operators("""
            x,            y,            z
            -x,            y,            z
            x,           -y,          1/2+z
            -x,           -y,          1/2+z
        """, offsets=[("+0", "+0", "+0"), ("+1/2", "+1/2", "+0")])
    elif spaceg == "C1c1" or spaceg == "9":
        symops = symmetry_operators("""
            x,            y,            z
            x,           -y,          1/2+z
        """, offsets=[("+0", "+0", "+0"), ("+1/2", "+1/2", "+0")])
    elif spaceg == "P1" or spaceg == "1":
        symops = symmetry_operators("""
            x,            y,            z
        """)
    elif spaceg == "Pbn2_1":
        # there are two #33 (?)
        symops = symmetry_operators("""
            x,            y,            z
            1/2-x,        1/2+y,          z
            1/2+x,        1/2-y,        1/2+z
            -x,           -y,          1/2+z
        """)
    elif spaceg == "Pna2_1":
        # there are two #33 (?)
        symops = symmetry_operators("""
            x,            y,            z
            1/2-x,        1/2+y,        1/2+z
            1/2+x,        1/2-y,          z
            -x,           -y,          1/2+z
        """)
    return symops

if __name__ == "__main__":
    for op in operations("Cmc2_1"):
        print(op(x=0.0, y=0.3, z=0.6))
