# coding: utf-8
desc = { "ref": { },
         "brief": "Rigid rotor (Quaternion).",
         "usage": "No options.",
         }

from logging import getLogger

import numpy as np
from scipy.spatial.transform import Rotation as R

import genice2.formats
from genice2.decorators import timeit, banner


class Format(genice2.formats.Format):
    """
The positions and orientations (in quaternions) of the water molecules are output in @NX4A format.
No options available.
    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)


    def hooks(self):
        return {5:self.Hook5}


    @timeit
    @banner
    def Hook5(self, ice):
        "Output water molecules as rigid rotors (Quaternion)."
        logger = getLogger()
        cellmat = ice.repcell.mat
        s = ""
        if cellmat[1,0] == 0 and cellmat[2,0] == 0 and cellmat[2,1] == 0:
            s += "@BOX3\n"
            s += "{0} {1} {2}\n".format(cellmat[0,0]*10,cellmat[1,1]*10,cellmat[2,2]*10)
        else:
            s += "@BOX9\n"
            for d in range(3):
                s += "{0} {1} {2}\n".format(cellmat[0,d]*10,cellmat[1,d]*10,cellmat[2,d]*10)
        s += "@NX4A\n"
        s += "{0}\n".format(len(ice.reppositions))
        for pos,rot in zip(ice.reppositions, ice.rotmatrices):
            position = np.dot(pos, cellmat)*10   #in Angstrom
            # quat     = rigid.rotmat2quat(rot.transpose())
            quat    = np.roll(R.from_matrix(rot).as_quat(),1) * np.array([-1,+1,-1,+1])
            s += " ".join([f"{v:9.4f}" for v in [*position, *quat]]) + "\n"
        s = "\n".join(ice.doc) + "\n" + s
        self.output = s
