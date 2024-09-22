# coding: utf-8
from genice2.decorators import timeit, banner
import genice2.formats
from logging import getLogger
from genice2 import rigid
import numpy as np
desc = {"ref": {},
        "brief": "Rigid rotor (Quaternion).",
        "usage": "No options.",
        }


class Format(genice2.formats.Format):
    """
The positions and orientations (in quaternions) of the water molecules are output in @NX4A format.
No options available.
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def hooks(self):
        return {5: self.Hook5}

    @timeit
    @banner
    def Hook5(self, ice):
        "Output water molecules as rigid rotors (Quaternion)."
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
        s += "@NX4A\n"
        s += "{0}\n".format(len(ice.reppositions))
        for pos, rot in zip(ice.reppositions, ice.rotmatrices):
            position = pos @ cellmat * 10  # in Angstrom
            quat = rigid.rotmat2quat(rot.transpose())
            s += "{0:9.4f} {1:9.4f} {2:9.4f}  {3:9.4f} {4:9.4f} {5:9.4f} {6:9.4f}\n".format(
                position[0], position[1], position[2], quat[0], quat[1], quat[2], quat[3])
        s = "\n".join(ice.doc) + "\n" + s
        self.output = s
