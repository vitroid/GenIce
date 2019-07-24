#!/usr/bin/python
"""
Crystallographic data of ice VIII

Kuhs, W. F., Finney, J. L., Vettier, C. & Bliss, D. V. Structure and hydrogen ordering in ices VI, VII, and VIII by neutron powder diffraction. J. Chem. Phys. 81, 3612–3623 (1998).
"""

atoms="""
O 0   0.25       0.1071(12)
D 0   0.4157(15) 0.1935(8)
"""

# space group: I4_1 /a m d No. 141

symops="""
     x,            y,            z
     -x,          1/2-y,          z
    3/4-y,        1/4+x,        3/4+z
    3/4+y,        3/4-x,        1/4+z
     -x,            y,            z
      x,          1/2-y,          z
    1/4+y,        1/4+x,        3/4+z
    1/4-y,        3/4-x,        1/4+z
 
     -x,           -y,           -z
      x,          1/2+y,         -z
    3/4+y,        1/4-x,        3/4-z
    3/4-y,        3/4+x,        1/4-z
      x,           -y,           -z
     -x,          1/2+y,         -z
    1/4-y,        1/4-x,        3/4-z
    1/4+y,        3/4+x,        1/4-z

""".translate({ord(','):''})

# add +1/2, +1/2, +1/2
lines = ""
for line in symops.split("\n"):
    cols = line.split()
    if len(cols) == 3:
        line = " ".join([x + "+1/2" for x in cols]) + "\n"
    lines += line

symops += lines


a=4.656 / 10.0 #nm
b=a
c=6.775 / 10.0 #nm
A=90
B=90
C=90

from genice.cell import cellvectors
cell  = cellvectors(a,b,c,A,B,C)

# helper routines to make from CIF-like data
from genice import CIF
atomd = CIF.atomdic(atoms)
sops  = CIF.symmetry_operators(symops)
waters, pairs = CIF.waters_and_pairs(cell, atomd, sops, rep=(2,2,2))

# All the pairs are fixed.
fixed = pairs

import numpy as np
density = 18*len(waters)/6.022e23 / (np.linalg.det(cell)*1e-21)

coord = "relative"
