# coding: utf-8

desc = { "ref": { },
         "brief": "Rigid rotor (Euler angle).",
         "usage": "No options.",
         }

from logging import getLogger

import numpy as np
from scipy.spatial.transform import Rotation as R

from genice2.decorators import timeit, banner
import genice2.formats


class Format(genice2.formats.Format):
    """
The positions and orientations (in Euler angles) of the water molecules are output in @NX3A format.
No options available.
    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)


    def hooks(self):
        return {5:self.Hook5}


    @timeit
    @banner
    def Hook5(self, ice):
        "Output water molecules as rigid rotors (Euler)."
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
        s += "@NX3A\n"
        s += "{0}\n".format(len(ice.reppositions))
        for pos,rot  in zip(ice.reppositions, ice.rotmatrices):
            position = np.dot(pos,cellmat)*10   #in Angstrom
            # euler = rigid.quat2euler(rigid.rotmat2quat(rot.transpose()))
            euler = R.from_matrix(rot.T).as_euler('ZXZ')
            euler = euler[[1,0,2]]
            euler[1:] -= np.floor(euler[1:]/(2*np.pi))*2*np.pi
            s += " ".join([f"{v:9.4f}" for v in [*position, *euler]]) + "\n"
        s = "\n".join(ice.doc) + "\n" + s
        self.output = s
