import logging

import numpy as np

from genice import rigid
from genice import yaplotlib as yp

def run(lattice, water_type="TIP3P", guests=[]):
    """
    Yaplot format.
    defined in https://github.com/vitroid/Yaplot
    """
    logger = logging.getLogger()
    lattice.stage1()   #replicate the unit cell
    lattice.stage2()   #prepare random graph

    logger.info("Output HBN in Yaplot format.")

    undir = lattice.graph.to_undirected()

    s = ""
    for i in lattice.graph.nodes_iter():
        pos = np.dot( lattice.reppositions[i], lattice.repcell )
        if i in lattice.dopants:
            s += yp.Color(6)
        elif 4 == len(undir.neighbors(i)):
            s += yp.Color(3)
        else:
            logger.debug("Z({1})={0}".format(undir.neighbors(i),i))
            s += yp.Color(5)
        s += yp.Layer(1)
        s += yp.Size(0.06)
        s += yp.Circle(pos)
        s += yp.Layer(2)
        if i in lattice.dopants:
            s += yp.Text(pos, "{0} ({1})".format(i, lattice.dopants[i]))
        else:
            s += yp.Text(pos, "{0}".format(i))
    s += yp.Color(2)
    for i,j in lattice.graph.edges_iter(data=False):
        s1 =lattice.reppositions[i]
        s2 =lattice.reppositions[j]
        d = s2-s1
        d -= np.floor( d + 0.5 )
        s2 = s1 + d
        s += yp.Layer(3)
        s += yp.Size(0.03)
        s += yp.Arrow(np.dot(s1,lattice.repcell),np.dot(s2,lattice.repcell))

    if not lattice.test2:
        # Failed to build the undirected graph obeying the ice rule.
        s = '#' + "\n#".join(lattice.doc) + "\n" + s
        print(s)
        return

    lattice.stage3()   #Make an ice graph
    lattice.stage4()   #Depolarize
    lattice.stage5()   #Orientation
    lattice.stage6(water_type)  #Water atoms
    lattice.stage7(guests)      #Guest atoms
    logger.info("Output water molecules in Yaplot format.")
    logger.info("Total number of atoms: {0}".format(len(lattice.atoms)))
    network = s
    s = lattice.yapresult
    s += yp.Layer(4)
    s += yp.ArrowType(1)
    H = []
    O  = ""
    for atom in lattice.atoms:
        resno, resname, atomname, position = atom
        if resno == 0:
            if O is not "":
                s += yp.Color(3)
                s += yp.Size(0.03)
                s += yp.Circle(O)
                if len(H):
                    s += yp.Line(O,H[0])
                    s += yp.Line(O,H[1])
                    s += yp.Size(0.01)
                    s += yp.Circle(H[0])
                    s += yp.Circle(H[1])
            H = []
            O = ""
        if "O" in atomname:
            O = position
        elif "H" in atomname:
            H.append(position)
        else:
            s += yp.Color(4)
            s += yp.Size(0.04)
            s += yp.Circle(position)
    if O is not "":
        s += yp.Circle(O)
        if len(H):
            s += yp.Line(O,H[0])
            s += yp.Line(O,H[1])
    s += network
    s = '#' + "\n#".join(lattice.doc) + "\n" + s
    print(s)
