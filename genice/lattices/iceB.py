#!/usr/bin/python

"""
**************************************************************************

BÁEZ, L. A. & CLANCY, P. Phase equilibria in extended simple point charge ice‐water systems. J. Chem. Phys. 103, 9744–9755 (1998).

**************************************************************************
"""

atoms="""
O1 0.2382 0.7981 0.3422
O2 0.5000 0.5000 0.0446
H1 -0.2603 0.9742 0.4613
H2 0.1446 0.8545 0.2013
H3 0.4108 0.6151 0.1577
"""

# P2_1 2_1 2
symops="""
     x            y            z
    1/2+x        1/2-y         -z
    1/2-x        1/2+y         -z
     -x           -y            z
"""

# in nm
a,b,c = 0.7176, 0.4428, 0.5040

from genice.cell import cellvectors
cell  = cellvectors(a,b,c)

# helper routines to make from CIF-like data
from genice import CIF
atomd = CIF.atomdic(atoms)
sops  = CIF.symmetry_operators(symops)
# the unit cell is too small to handle; multiply (2,2,2)
waters, fixed = CIF.waters_and_pairs(cell, atomd, sops, rep=(2,2,2))

# set pairs in this way for hydrogen-ordered ices.
pairs = fixed

import numpy as np
density = 18*len(waters)/6.022e23 / (np.linalg.det(cell)*1e-21)

coord = "relative"
