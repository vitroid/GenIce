# coding: utf-8

desc={"ref": {},
      "brief": "MDView file (in Angdtrom).",
      "usage": "No options available."
      }

if __name__[-6:] == 'mdv_au':
    desc["brief"] = "MDView file (in au)."

import numpy as np
from logging import getLogger

au = 0.052917721092 # nm
AA = 0.1 # nm

import genice.formats
class Format(genice.formats.Format):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)


    def hooks(self):
        return {7:self.hook7}


    def hook7(self, ice):
        logger = getLogger()
        logger.info("Hook7: Output in MDView format.")
        logger.info("  Total number of atoms: {0}".format(len(ice.atoms)))
        if __name__[-6:] == 'mdv_au':
            conv = 1.0 / au
        else:
            conv = 10.0
        cellmat = ice.repcell.mat
        s = ""
        s += "# {0} {1} {2}\n".format(cellmat[0,0]*conv, cellmat[1,1]*conv, cellmat[2,2]*conv)
        s += "-center 0 0 0\n"
        s += "-fold\n"
        s += "{0}\n".format(len(ice.atoms))
        for atom in ice.atoms:
            molorder, resname, atomname, position, order = atom
            s += "{0:5} {1:9.4f} {2:9.4f} {3:9.4f}\n".format(atomname,*(position[:3]*conv))
        self.output = s
        logger.info("Hook7: end.")
