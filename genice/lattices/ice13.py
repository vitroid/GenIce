#!/usr/bin/python

# Crystallographic data of ice XIII
# Salzmann, C. G., Radaelli, P., Hallbrucker, A. & Mayer, E. The Preparation and Structures of Hydrogen Ordered Phases of Ice. Science 311, 1758–1761 (2006).

atoms="""
O1 0.2541(6)  0.5629(5) 0.2517(5)
O2 0.4771(6)  0.7992(5) 0.4089(5)
O3 0.0503(6)  0.8082(6) 0.0941(5)
O4 0.2613(5)  0.4045(6) 0.4992(5)
O5 0.2113(4)  0.4029(5) 0.0034(5)
O6 0.4147(5)  0.1103(7) 0.2336(4)
O7 0.1245(5)  0.1142(6) 0.2643(4)
D8 0.3444(4)  0.6427(5) 0.3008(3)
D10 0.2458(5) 0.4942(5) 0.3299(5)
D13 0.1074(4) 0.7187(5) 0.1563(4)
D16 0.4820(4) 0.9075(5) 0.3558(4)
D18 0.5763(5) 0.7499(5) 0.4437(4)
D19 0.9486(5) 0.7508(5) 0.0478(4)
D21 0.2372(3) 0.4543(5) 0.0989(4)
D24 0.3043(4) 0.4904(6) 0.5777(4)
D26 0.1708(4) 0.3555(6) 0.5137(4)
D27 0.3072(4) 0.3737(6) 0.9904(3)
D29 0.0781(4) 0.0194(6) 0.1989(4)
D30 0.3250(5) 0.1374(5) 0.2554(5)
D32 0.3823(5) 0.0496(6) 0.1467(5)
D35 0.0509(4) 0.2082(6) 0.2548(5)
"""

# space group: P2_1/a No. 14
# http://img.chem.ucl.ac.uk/sgp/large/014dy1.htm

symops="""
    x            y            z
    1/2-x        1/2+y         -z
 
     -x           -y           -z
    1/2+x        1/2-y          z
"""

a=9.2417 / 10.0 #nm
b=7.4724 / 10.0 #nm
c=10.2970 / 10.0 #nm
A=90
B=109.6873
C=90

from genice.cell import cellvectors
cell  = cellvectors(a,b,c,A,B,C)

# helper routines to make from CIF-like data
from genice import CIF
atomd = CIF.atomdic(atoms)
sops  = CIF.symmetry_operators(symops)
waters, fixed = CIF.waters_and_pairs(cell, atomd, sops)

# set pairs in this way for hydrogen-ordered ices.
pairs = fixed

import numpy as np
density = 18*len(waters)/6.022e23 / (np.linalg.det(cell)*1e-21)

coord = "relative"
