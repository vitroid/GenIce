# coding: utf-8
desc={"ref": {"exyz": "https://open-babel.readthedocs.io/en/latest/FileFormats/Extended_XYZ_cartesian_coordinates_format.html"},
      "brief": "Extended XYZ format.",
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


    def hook7(self, lattice):
        logger = getLogger()
        logger.info("Hook7: Output in extended XYZ format.")
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
        logger.info("Hook7: end.")
