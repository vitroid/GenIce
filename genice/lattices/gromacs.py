# -*- coding:utf-8 -*-

"""
Load a gromacs file as a lattice.

Usage: genice gromacs[filename]            Regards Ow and Hw to be of a water molecule.
       genice gromacs[filename:Ow]         Specify the name of O atom (and ignore H positions)
       genice gromacs[filename:Ow:Hw[12]]  Specify the names of O and H atoms. (Regexp is accepted.)

"""


import numpy as np
import re
import logging
import genice.pairlist as pl

def argparser(arg):
    global waters, cell, celltype, coord, density, pairs
    logger = logging.getLogger()
    args = arg.split(":")
    assert 0 < len(args) <= 3, __doc__
    if len(args) == 1:
        O = "Ow"
        H = "Hw"
    elif len(args) == 2:
        O = args[1]
        H = None
    else:
        O = args[1]
        H = args[2]
    filename = args[0]

    file = open(filename)
    file.readline()
    natom = int(file.readline())
    hatoms = []
    waters = []
    for i in range(natom):
        line = file.readline()
        # resid = int(line[0:5])
        # resna = line[5:10]
        atomname = line[10:15].replace(' ', '')
        # atomid = int(line[15:20])
        pos = np.array([float(x) for x in line[20:].split()[:3]]) #drop velocity
        if atomname == O:
            waters.append(pos)
        elif H is not None and re.fullmatch(H, atomname):
            hatoms.append(pos)
        else:
            logger.info("Skip {0}".format(atomname))
    c = [float(x) for x in file.readline().split()]
    if len(c) == 3:
        cell = np.array([[c[0],0.,0.],
                         [0.,c[1],0.],
                         [0.,0.,c[2]]])
    else:
        cell = np.array([[c[0],c[3],c[4]],
                         [c[5],c[1],c[6]],
                         [c[7],c[8],c[2]]])
    celltype = 'triclinic'
    coord = 'absolute'
    density = len(waters) / (np.linalg.det(cell)*1e-21) * 18 / 6.022e23

    if len(hatoms) > 0:
        celli = np.linalg.inv(cell)
        # relative coord
        rh = [np.dot(x, celli) for x in hatoms]
        ro = [np.dot(x, celli) for x in waters]
        grid = pl.determine_grid(cell, 0.245)
        pairs0 = pl.pairlist_fine_hetero(ro, rh, 0.245, cell, grid, distance=False)
        # remove intramolecular OHs
        pairs = []
        for o,h in pairs0:
            if h == o*2 or h == o*2+1:
                # adjust oxygen positions
                dh = rh[h] - ro[o]
                dh -= np.floor(dh + 0.5)
                waters[o] += np.dot(dh, cell)*1./16.
            else:
                # register a new intermolecular pair
                pairs.append((h//2, o))
        logger.debug("# of pairs: {0} {1}".format(len(pairs),len(waters)))
        
    
    

# default. Do nothing (leading to an error).

