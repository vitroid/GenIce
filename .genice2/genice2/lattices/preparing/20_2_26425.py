from genice2.cell import cellvectors
import genice2.lattices
import numpy as np
desc = {
    "ref": {
        "20_2_26425": "Engel 2018",
        "engel08": "Engel 2018"
    },
    "usage": "No options available.",
    "brief": "Hypothetical zeolitic ice",
    "test": ({"args": "",
             "options": "-r 2 2 2"}, )
}


class Lattice(genice2.lattices.Lattice):
    def __init__(self):
        self.cell = np.array([
            [7.356626, 2.559948, 0.030551],
            [4.253264, 6.305892, 0.011851],
            [-0.016176, 0.001573, 3.176074],
        ])
        self.waters = np.array([
            [-0.139504, 0.427152, -0.089349],
            [0.056636, -0.446229, 0.412664],
            [-0.441797, 0.066077, -0.211429],
            [0.359048, -0.085215, 0.285567],
            [0.195383, -0.246391, 0.592419],
            [-0.269726, 0.219171, 0.091624],
        ])
        self.coord = 'relative'
