# coding: utf-8
import numpy as np
#United-atom THF model with a dummy site

from logging import getLogger
import genice.molecules

class Molecule(genice.molecules.Molecule):
    def __init__(self):
        self.sites = np.array([[ 0.,       -0.119625,  0.      ],
                               [ 0.116284, -0.039705,  0.      ],
                               [ 0.076453,  0.107915,  0.      ],
                               [-0.076453,  0.107915,  0.      ],
                               [-0.116284, -0.039705,  0.      ],
                               [ 0.0,       0.0,       0.      ] ])

        self.atoms = ["O","C","C","C","C","."]
        self.labels = ["O","CA","CB","CB","CA","CM"]
        self.name = "THF"
