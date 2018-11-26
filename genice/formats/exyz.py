# coding: utf-8
"""
extended XMol file format (.exyz)

http://open-babel.readthedocs.io/en/latest/FileFormats/Extended_XYZ_cartesian_coordinates_format.html
"""

import numpy as np

from genice import rigid

def hook7(lattice):
    lattice.logger.info("Hook7: Output in extended XYZ format.")
    s = ""
    s += "{0}\n".format(len(lattice.atoms))
    s += "%PBC\n"
    for atom in lattice.atoms:
        molorder, resname, atomname, position, order = atom
        s += "{0:>4}{1:15.5f}{2:15.5f}{3:15.5f}\n".format(atomname,position[0]*10,position[1]*10,position[2]*10)
    s = '#' + "\n#".join(lattice.doc) + "\n" + s
    s += "\n"
    s += "Vector1 {0:15.5f}{1:15.5f}{2:15.5f}\n".format(*lattice.repcell.mat[0]*10)
    s += "Vector2 {0:15.5f}{1:15.5f}{2:15.5f}\n".format(*lattice.repcell.mat[1]*10)
    s += "Vector3 {0:15.5f}{1:15.5f}{2:15.5f}\n".format(*lattice.repcell.mat[2]*10)
    s += "Offset  {0:15.5f}{1:15.5f}{2:15.5f}\n".format(0.,0.,0.)
    print(s,end="")
    lattice.logger.info("Hook7: end.")


hooks = {7:hook7}
