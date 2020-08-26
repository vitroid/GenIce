# coding: utf-8

import numpy as np

class Molecule():
    """
    Base class of molecules
    """
    def __init__(self, **kwargs):
        assert len(kwargs) == 0

        # sites: positions of interaction sites relative to a center of molecule
        # a numpy array of rank (N, 3) where N is number of sites
        self.sites = np.zeros([1,3])

        # Labels of the interaction sites.
        self.labels = ["Me",]

        # the name that represents the molecule. It is necessary for Gromacs format.
        self.name = "MET"

        # Define a list of atom names when the molecule contains virtual sites.
        self.atoms = None
