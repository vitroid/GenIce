import logging

import numpy as np

def run(lattice, water_type="TIP3P", guests=[]):
    """
    re-make python module for GenIce
    """
    logger = logging.getLogger()
    lattice.stage1()   #replicate the unit cell
    lattice.stage2()   #prepare random graph
    lattice.stage3()   #Make an ice graph
    lattice.stage4()   #Depolarize
    lattice.stage5()   #Orientation
    lattice.stage6(water_type)  #Water atoms
    lattice.stage7(guests)      #Guest atoms
    logger.info("Total number of atoms: {0}".format(len(lattice.atoms)))
    logger.info("Output as a python module.")
    s = ""
    s += '"""\n'
    s += "\n".join(lattice.doc) + "\n"
    s += '"""\n'
    s += "bondlen={0}\n".format(lattice.bondlen)
    s += "coord='relative'\n"
    if lattice.repcell[1,0] == 0 and lattice.repcell[2,0] == 0 and lattice.repcell[2,1] == 0:
        s += "celltype='rect'\n"
        s += "cell='{0} {1} {2}'\n".format(lattice.repcell[0,0],lattice.repcell[1,1],lattice.repcell[2,2])
    else:
        s += "celltype='triclinic'\n"
        s += "cell='{0} {1} {2} {3} {4} {5} {6} {7} {8}'\n".format(*lattice.repcell[0],*lattice.repcell[1], *lattice.repcell[2])
    s += "density={0}\n".format(lattice.density)
    s += "waters=\"\"\"\n"
    for i in range(len(lattice.reppositions)):
        position = lattice.reppositions[i]
        s += "{0:9.4f} {1:9.4f} {2:9.4f}\n".format(position[0],position[1],position[2])
    s += "\"\"\"\n\n"
    print(s,end="")
