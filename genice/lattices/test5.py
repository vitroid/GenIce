import numpy as np


#Inherit HS1
from genice.lattices.HS1 import *


#This is called at the end of Lattice.__init__()
def dopeIonsToUnitCell(lat):
    """
    Special post-process for inherited module
    """
    #NH4F
    #3 --> NH4
    #39 --> F
    #it is a way to put ions at the fixed positions in the unit cell.
    #They will be replicated when the cell is replicated.
    lat.logger.info("  Start doping.")
    lat.graph.cationize(3)
    lat.graph.cationize(63)
    lat.graph.anionize(39)
    lat.graph.anionize(19)
    lat.dopants = {3:"Nc", 63:"Nc", 39:"Br", 19:"Br"}

    #Show info on the environment
    #for the case a large molecular ions are doped.
    lat.dopants_info()

    #add groups
    lat.add_group(cage=9, group="ethyl", root=3)
    lat.add_group(cage=11, group="propyl", root=3)
    lat.add_group(cage=13, group="butyl", root=3)
    lat.add_group(cage=7, group="pentyl", root=3)
    lat.add_group(cage=8, group="3-methylbutyl", root=63)
    lat.add_group(cage=10, group="3,3-dimethylbutyl", root=63)
    lat.add_group(cage=12, group="2,3-dimethylbutyl", root=63)
    lat.add_group(cage=6, group="2,2-dimethylpropyl", root=63)
                  
    lat.logger.info("  end.")
