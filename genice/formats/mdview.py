# coding: utf-8

desc={"ref": {},
      "brief": "MDView file (in Angdtrom).",
      "usage": "No options available."
      }

if __name__[-6:] == 'mdv_au':
    desc["brief"] = "MDView file (in au)."

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
    cellmat = lattice.repcell.mat
    s = ""
    s += "# {0} {1} {2}\n".format(cellmat[0,0]*conv, cellmat[1,1]*conv, cellmat[2,2]*conv)
    s += "-center 0 0 0\n"
    s += "-fold\n"
    s += "{0}\n".format(len(lattice.atoms))
    for atom in lattice.atoms:
        molorder, resname, atomname, position, order = atom
        s += "{0:5} {1:9.4f} {2:9.4f} {3:9.4f}\n".format(atomname,*(position[:3]*conv))
    print(s,end="")
    lattice.logger.info("Hook7: end.")


hooks = {7:hook7}
