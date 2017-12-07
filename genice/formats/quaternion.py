# coding: utf-8
"""
Rigid water molecule 
"""

import numpy as np

from genice import rigid


def hook5(lattice):
    lattice.logger.info("Hook5: Output water molecules as rigid rotors (Quaternion).")
    cellmat = lattice.repcell.mat
    s = ""
    if cellmat[1,0] == 0 and cellmat[2,0] == 0 and cellmat[2,1] == 0:
        s += "@BOX3\n"
        s += "{0} {1} {2}\n".format(cellmat[0,0]*10,cellmat[1,1]*10,cellmat[2,2]*10)
    else:
        s += "@BOX9\n"
        for d in range(3):
            s += "{0} {1} {2}\n".format(cellmat[0,d]*10,cellmat[1,d]*10,cellmat[2,d]*10)
    s += "@NX4A\n"
    s += "{0}\n".format(len(lattice.reppositions))
    for i in range(len(lattice.reppositions)):
        position = np.dot(lattice.reppositions[i],cellmat)*10   #in Angstrom
        quat     = rigid.rotmat2quat(lattice.rotmatrices[i].transpose())
        s += "{0:9.4f} {1:9.4f} {2:9.4f}  {3:9.4f} {4:9.4f} {5:9.4f} {6:9.4f}\n".format(position[0],
                                                                        position[1],
                                                                        position[2],
                                                                        quat[0],
                                                                        quat[1],
                                                                        quat[2],
                                                                        quat[3])
    s = "\n".join(lattice.doc) + "\n" + s
    print(s,end="")
    lattice.logger.info("Hook5: end.")


hooks = {5:hook5}
