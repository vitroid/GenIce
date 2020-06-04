#!/usr/bin/python
"""
Usage: genice ice1hte
"""
from logging import getLogger
import numpy as np

def usage():
    logger = getLogger()
    logger.info(__doc__)

desc={"ref": {"Ih": "Page 11 of the Supplemenrary Material of P. Teeratchanan and A. Hermann, Computational phase diagrams of noble gas hydrates under pressure, J. Chem. Phys. 143, 154507 (2015); https://doi.org/10.1063/1.4933371",
},
      "usage": usage(),
      "brief": "Filled ice Ih by Teeratchanan (Hydrogen disordered). (Positions of guests are supplied.)"
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
    global pairs, fixed, waters, coord, density, cell, cagetype, cagepos
    logger = getLogger()

    # Ref. Ih
    atoms="""
O1 0.0000 0.6699 0.0488 
O2 0.0000 0.3377 -0.0557
Ne1 0.0000 0.0013 0.7539
    """

    # Ref. Ih
    # space group: Cmc2_1
    symops="""
      x,            y,            z
     -x,            y,            z
      x,           -y,          1/2+z
     -x,           -y,          1/2+z
      x+1/2,            y+1/2,            z
     -x+1/2,            y+1/2,            z
      x+1/2,           -y+1/2,          1/2+z
     -x+1/2,           -y+1/2,          1/2+z
    """.replace(',', ' ')

    # Ref. Ih
    a=4.568 / 10.0 #nm
    b=7.980 / 10.0 #nm
    c=6.894 / 10.0 #nm

    from genice.cell import cellvectors
    cell  = cellvectors(a,b,c)

    # helper routines to make from CIF-like data
    from genice import CIF
    atomd = CIF.atomdic(atoms)
    atoms = CIF.fullatoms(atomd, CIF.symmetry_operators(symops))

    cagetype = []
    cagepos  = []
    for name, pos in pick_atoms(atoms, ("Ne1",), repeat=(2,1,1)):
        cagetype.append(name)
        cagepos.append(pos)
    
    waters, pairs = CIF.waters_and_pairs(cell, atomd, CIF.symmetry_operators(symops), rep=(2,1,1))

    density = 18*len(waters)/6.022e23 / (np.linalg.det(cell)*1e-21)

    coord = "relative"
