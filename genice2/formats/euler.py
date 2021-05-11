# coding: utf-8

import genice2.formats
from genice2.decorators import timeit, banner
from logging import getLogger
from genice2 import rigid
import numpy as np
desc = {"ref": {},
        "brief": "Rigid rotor (Euler angle).",
        "usage": "No options.",
        }


class Format(genice2.formats.Format):
    """
The positions and orientations (in Euler angles) of the water molecules are output in @NX3A format.
No options available.
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def hooks(self):
        return {5: self.Hook5}

    @timeit
    @banner
    def Hook5(self, ice):
        "Output water molecules as rigid rotors (Euler)."
        logger = getLogger()
        cellmat = ice.repcell.mat
        s = ""
        if cellmat[1, 0] == 0 and cellmat[2, 0] == 0 and cellmat[2, 1] == 0:
            s += "@BOX3\n"
            s += "{0} {1} {2}\n".format(cellmat[0, 0] *
                                        10, cellmat[1, 1] * 10, cellmat[2, 2] * 10)
        else:
            s += "@BOX9\n"
            for d in range(3):
                s += "{0} {1} {2}\n".format(cellmat[0, d] *
                                            10, cellmat[1, d] * 10, cellmat[2, d] * 10)
        s += "@NX3A\n"
        s += "{0}\n".format(len(ice.reppositions))
        for pos, rot in zip(ice.reppositions, ice.rotmatrices):
            position = pos @ cellmat * 10  # in Angstrom
            euler = rigid.quat2euler(rigid.rotmat2quat(rot.transpose()))
            s += "{0:9.4f} {1:9.4f} {2:9.4f}  {3:9.4f} {4:9.4f} {5:9.4f}\n".format(
                position[0], position[1], position[2], euler[0], euler[1], euler[2])
        s = "\n".join(ice.doc) + "\n" + s
        self.output = s
