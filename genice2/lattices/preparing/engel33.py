from genice2.cell import cellvectors
import genice2.lattices
import numpy as np
desc = {
    "ref": {
        "14_2_48453": "Engel 2018",
        "engel33": "Engel 2018"
    },
    "usage": "No options available.",
    "brief": "Hypothetical zeolitic ice",
    "test": ({"args": "",
             "options": "-r 2 2 2"},)
}


class Lattice(genice2.lattices.Lattice):
    def __init__(self):
        self.cell = np.array([
            [4.274028, -0.495899, -0.007409],
            [0.69394, 6.104765, 0.009168],
            [-0.074983, 0.005107, 6.682336],
        ])
        self.waters = np.array([
            [0.447519, 0.104885, -0.25555],
            [0.434958, -0.054775, 0.184694],
            [-0.053005, 0.445498, 0.24353],
            [-0.07881, -0.395455, -0.316209],
            [0.449664, 0.2932, 0.396693],
            [0.417874, -0.241958, -0.467209],
            [-0.048131, 0.257569, -0.104188],
            [-0.069137, -0.206786, 0.031111],
        ])
        self.coord = 'relative'
