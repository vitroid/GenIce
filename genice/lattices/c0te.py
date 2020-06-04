#!/usr/bin/python
"""
Usage: genice c0te
"""

from logging import getLogger
import numpy as np

def usage():
    logger = getLogger()
    logger.info(__doc__)

desc={"ref": {"C0": "Page 11 of the Supplemenrary Material of P. Teeratchanan and A. Hermann, Computational phase diagrams of noble gas hydrates under pressure, J. Chem. Phys. 143, 154507 (2015); https://doi.org/10.1063/1.4933371",
},
      "usage": usage(),
      "brief": "Filled ice C0 by Teeratchanan (Hydrogen-disordered.) (Positions of guests are supplied.)"
      }



def argparser(arg):
    global pairs, fixed, waters, coord, density, cell, cagetype, cagepos
    logger = getLogger()

    # Ref. 2atom
    atoms="""
    O1 0.2342 0.4721 0.8019
    O2 0.7648 0.5306 0.2941
    Ne1 -0.0647 0.7868 0.7669
    """

    # Ref. 
    # space group: P3_2
    symops="""
      x,            y,            z
-y,x-y,z+2/3
-x+y,-x,z+1/3
    """.replace(',', ' ')

    # Ref. 2cell
    a=6.177 / 10.0 #nm
    c=6.054 / 10.0 #nm
    C=120.0

    from genice.cell import cellvectors
    cell  = cellvectors(a,a,c,C=C)

    # helper routines to make from CIF-like data
    from genice import CIF
    atomd = CIF.atomdic(atoms)
    atoms = CIF.fullatoms(atomd, CIF.symmetry_operators(symops))

    cagetype = []
    cagepos  = []
    for atomname, pos in atoms:
        if atomname == "Ne1":
            cagetype.append(atomname)
            cagepos.append(pos)
            
    waters, pairs = CIF.waters_and_pairs(cell, atomd, CIF.symmetry_operators(symops))

    density = 18*len(waters)/6.022e23 / (np.linalg.det(cell)*1e-21)

    coord = "relative"
