# coding: utf-8
from genice2.molecules import serialize
import genice2.formats
from genice2.decorators import timeit, banner
from genice2 import rigid
from logging import getLogger
import numpy as np
desc = {
    "ref": {
        "exyz": "https://open-babel.readthedocs.io/en/latest/FileFormats/Extended_XYZ_cartesian_coordinates_format.html"},
    "brief": "Extended XYZ format.",
    "usage": "No options available."}


class Format(genice2.formats.Format):
    """
The atomic positions of the molecules are output in an extended XYZ format.
No options available.
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def hooks(self):
        return {7: self.Hook7}

    @timeit
    @banner
    def Hook7(self, ice):
        "Output in extended XYZ format."
        logger = getLogger()
        atoms = []
        for mols in ice.universe:
            atoms += serialize(mols)
        s = ""
        s += "{0}\n".format(len(atoms))
        s += "%PBC\n"
        for atom in atoms:
            molorder, resname, atomname, position, order = atom
            s += "{0:>4}{1:15.5f}{2:15.5f}{3:15.5f}\n".format(
                atomname, position[0] * 10, position[1] * 10, position[2] * 10)
        s = '#' + "\n#".join(ice.doc) + "\n" + s
        s += "\n"
        s += "Vector1 {0:15.5f}{1:15.5f}{2:15.5f}\n".format(
            *ice.repcell.mat[0] * 10)
        s += "Vector2 {0:15.5f}{1:15.5f}{2:15.5f}\n".format(
            *ice.repcell.mat[1] * 10)
        s += "Vector3 {0:15.5f}{1:15.5f}{2:15.5f}\n".format(
            *ice.repcell.mat[2] * 10)
        s += "Offset  {0:15.5f}{1:15.5f}{2:15.5f}\n".format(0., 0., 0.)
        self.output = s
