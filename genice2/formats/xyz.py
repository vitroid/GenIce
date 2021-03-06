# coding: utf-8
"""
XMol file format (.xyz)
"""

import numpy as np

from genice2 import rigid
from logging import getLogger
from genice2.decorators import timeit, banner
import genice2.formats


class Format(genice2.formats.Format):
    """
The atomic positions of the molecules are output in a crude XYZ format.
No options available.
    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)


    def hooks(self):
        return {7:self.Hook7}


    @timeit
    @banner
    def Hook7(self, ice):
        "Output in XYZ format."
        logger = getLogger()
        s = ""
        s += "{0}\n".format(len(ice.atoms))
        s += "Generated by GenIce. https://github.com/vitroid/GenIce\n"
        for atom in ice.atoms:
            molorder, resname, atomname, position, order = atom
            s += "{0:5} {1:9.4f} {2:9.4f} {3:9.4f}\n".format(atomname,position[0]*10,position[1]*10,position[2]*10)
        s = '#' + "\n#".join(ice.doc) + "\n" + s
        self.output = s
