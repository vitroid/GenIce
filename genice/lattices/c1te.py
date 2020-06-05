#!/usr/bin/python
"""
Usage: genice c1te
"""
from logging import getLogger
import numpy as np

def usage():
    logger = getLogger()
    logger.info(__doc__)

desc={"ref": {"C1": "Page 12 of the Supplemenrary Material of P. Teeratchanan and A. Hermann, Computational phase diagrams of noble gas hydrates under pressure, J. Chem. Phys. 143, 154507 (2015); https://doi.org/10.1063/1.4933371"},
      "usage": usage(),
      "brief": "Hydrogen-ordered hydrogen hydrate C1 by Teeratchanan. (Positions of guests are supplied.)"
      }

def pick_atoms(atoms, names, repeat=(1,1,1)):
    nrep = np.array(repeat)
    for atomname, fracpos in atoms:
        if atomname in names:
            for x in range(repeat[0]):
                for y in range(repeat[1]):
                    for z in range(repeat[2]):
                        yield atomname, (fracpos+np.array([x,y,z]))/nrep


def argparser(arg):
    global pairs, fixed, waters, coord, density, cell, cagepos, cagetype
    logger = getLogger()

    # Ref. C1
    atoms="""
O1 0.2228 0.1966 0.0454
O2 0.5236 0.8974  0.1466
H1 0.8184 0.5344  0.3198
H2 0.2244 0.2188 0.2045
H3 0.4384 0.8815  0.1496
H4 0.5729 -0.0253 0.2297
Ne1 0.0000 0.0000 0.7361
    """

    # Ref. C1
    # space group: R-3 No. 148
    # in a rhombus cell
    symops="""
     x,            y,            z
     -y,           x-y,           z
    -x+y,          -x,            z
 
     -x,           -y,           -z
      y,          -x+y,          -z
     x-y,           x,           -z

      x+2/3,            y+1/3,           z+1/3
     -y+2/3,          x-y+1/3,           z+1/3
   -x+y+2/3,           -x+1/3,           z+1/3
 
     -x+2/3,           -y+1/3,          -z+1/3
      y+2/3,         -x+y+1/3,          -z+1/3
    x-y+2/3,            x+1/3,          -z+1/3

      x+1/3,            y+2/3,           z+2/3
     -y+1/3,          x-y+2/3,           z+2/3
   -x+y+1/3,           -x+2/3,           z+2/3
 
     -x+1/3,           -y+2/3,          -z+2/3
      y+1/3,         -x+y+2/3,          -z+2/3
    x-y+1/3,            x+2/3,          -z+2/3
    """.replace(',', ' ')

    # Ref. C1
    a=12.673 / 10.0 #nm
    c= 6.017 / 10.0 #nm
    C= 120

    from genice.cell import cellvectors
    cell  = cellvectors(a,a,c,C=C)

    # helper routines to make from CIF-like data
    from genice import CIF
    atomd = CIF.atomdic(atoms)
    atoms = CIF.fullatoms(atomd, CIF.symmetry_operators(symops))

    cagetype = []
    cagepos  = []
    for name, pos in pick_atoms(atoms, ("Ne1",)):
        cagetype.append(name)
        cagepos.append(pos)
        
    sops  = CIF.symmetry_operators(symops)
    waters, fixed = CIF.waters_and_pairs(cell, atomd, sops)

    # set pairs in this way for hydrogen-ordered ices.
    pairs = fixed

    density = 18*len(waters)/6.022e23 / (np.linalg.det(cell)*1e-21)

    coord = "relative"
