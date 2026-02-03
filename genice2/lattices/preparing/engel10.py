from genice2.cell import cellvectors
import genice2.lattices
import numpy as np
desc = {
    "ref": {
        "84_2_1419": "Engel 2018",
        "engel10": "Engel 2018"
    },
    "usage": "No options available.",
    "brief": "Hypothetical zeolitic ice"
}


class Lattice(genice2.lattices.Lattice):
    def __init__(self):
        self.cell = np.array([
            [4.212594, -0.029181, 0.000138],
            [0.008311, 4.213037, -4.5e-05],
            [0.000235, -7.1e-05, 7.499247],
        ])
        self.waters = np.array([
            [0.199671, 0.233915, -0.499994],
            [0.297268, -0.172515, 1e-06],
            [-0.149621, -0.255634, -0.500009],
            [-0.1882, 0.177552, -1e-06],
            [-0.46738, 0.48519, 0.2491],
            [-0.467389, 0.485201, -0.2491],
        ])
        self.coord = 'relative'
