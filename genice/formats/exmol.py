# coding: utf-8
desc={"ref": {"exmol": "http://libatoms.github.io/QUIP/io.html#extendedxyz"},
      "brief": "Extended XMol file format.",
      "usage": "No options available."
      }


import numpy as np
from logging import getLogger
from genice import rigid

import genice.formats
class Format(genice.formats.Format):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)


    def hooks(self):
        return {7:self.hook7}


    def hook7(self, ice):
        logger = getLogger()
        logger.info("Hook7: Output in extended XMol file format.")
        s = ""
        s += "{0}\n".format(len(ice.atoms))
        s += 'Lattice="{0:.3f} {1:.3f} {2:.3f} {3:.3f} {4:.3f} {5:.3f} {6:.3f} {7:.3f} {8:.3f}"\n'.format(*ice.repcell.mat.reshape((9,))*10)
        for atom in ice.atoms:
            molorder, resname, atomname, position, order = atom
            s += "{0:>4}{1:15.5f}{2:15.5f}{3:15.5f}\n".format(atomname,position[0]*10,position[1]*10,position[2]*10)
        s = '#' + "\n#".join(ice.doc) + "\n" + s
        print(s,end="")
        logger.info("Hook7: end.")
