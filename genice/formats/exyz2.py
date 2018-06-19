# coding: utf-8
"""
extended XMol file format (.exyz)

http://libatoms.github.io/QUIP/io.html#extendedxyz
"""

import numpy as np

from genice import rigid

def hook7(lattice):
    lattice.logger.info("Hook7: Output in extended XYZ format.")
    s = ""
    s += "{0}\n".format(len(lattice.atoms))
    s += 'Lattice="{0:.3f} {1:.3f} {2:.3f} {3:.3f} {4:.3f} {5:.3f} {6:.3f} {7:.3f} {8:.3f}"\n'.format(*lattice.repcell.mat.reshape((9,))*10)
    for atom in lattice.atoms:
        molorder, resname, atomname, position, order = atom
        s += "{0:>4}{1:15.5f}{2:15.5f}{3:15.5f}\n".format(atomname,position[0]*10,position[1]*10,position[2]*10)
    s = '#' + "\n#".join(lattice.doc) + "\n" + s
    print(s,end="")
    lattice.logger.info("Hook7: end.")


hooks = {7:hook7}
