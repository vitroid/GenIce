# coding: utf-8
"""
MDView file format
"""

import numpy as np

au = 0.052917721092 # nm
AA = 0.1 # nm

def hook7(lattice):
    lattice.logger.info("Hook7: Output in MDView format.")
    lattice.logger.info("  Total number of atoms: {0}".format(len(lattice.atoms)))
    if __name__[-6:] == 'mdv_au':
        conv = 1.0 / au
    else:
        conv = 10.0
    s = ""
    s += "-center 0 0 0\n"
    s += "-fold\n"
    s += "{0}\n".format(len(lattice.atoms))
    for i in range(len(lattice.atoms)):
        molorder, resname, atomname, position, order = lattice.atoms[i]
        s += "{0:5} {1:9.4f} {2:9.4f} {3:9.4f}\n".format(atomname,*(position[:3]*conv))
    s = '#' + "\n#".join(lattice.doc) + "\n" + s
    print(s,end="")
    lattice.logger.info("Hook7: end.")


hooks = {7:hook7}
