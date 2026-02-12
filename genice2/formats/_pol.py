
# coding: utf-8

from genice2.decorators import timeit, banner
import genice2.formats
from logging import getLogger
import numpy as np
desc = {
    "ref": {},
    "brief": "Polarization check.",
    "usage": "Polarization check.",
}


# It should be expressed as a function of distance.

class Format(genice2.formats.Format):
    """
Calculate the polarization.
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def hooks(self):
        return {4: self.Hook4}

    @timeit
    @banner
    def Hook4(self, ice):
        "Polarization check."
        
        logger = getLogger()
        
        polarization = 0
        for i, j in ice.digraph.edges():
            pos_i = ice.reppositions[i]
            pos_j = ice.reppositions[j]
            dvec = pos_j - pos_i
            dvec -= np.floor(dvec + 0.5)
            polarization += dvec

        self.output = f"{polarization[0]:.5f} {polarization[1]:.5f} {polarization[2]:.5f} polarization\n"
