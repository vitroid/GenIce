# coding: utf-8
"""
PovRay format.
"""

import numpy as np

from genice import rigid
from genice import povraylib as po


def hook6(lattice):
    lattice.logger.info("Output water molecules in PovRay format.")
    lattice.logger.info("Total number of atoms: {0}".format(len(lattice.atoms)))

    #prepare the 

        # network = self.s
    s = ""
    H = []
    O  = ""
    for atom in lattice.atoms:
        resno, resname, atomname, position = atom
        if resno == 0:
            if O is not "":
                s += po.Color(3)
                s += po.Size(0.03)
                s += po.Circle(O)
                if len(H):
                    s += po.Line(O,H[0])
                    s += po.Line(O,H[1])
                    s += po.Size(0.01)
                    s += po.Circle(H[0])
                    s += po.Circle(H[1])
            H = []
            O = ""
        if "O" in atomname:
            O = position
        elif "H" in atomname:
            H.append(position)
        else:
            s += po.Color(4)
            s += po.Size(0.04)
            s += po.Circle(position)
    if O is not "":
        s += po.Circle(O)
        if len(H):
            s += po.Line(O,H[0])
            s += po.Line(O,H[1])
    for i,j in lattice.graph.edges_iter(data=False):

    s += network
    s = '#' + "\n#".join(lattice.doc) + "\n" + s
    print(s)


hooks = {6:hook6}
