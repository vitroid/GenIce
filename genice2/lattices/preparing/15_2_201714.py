from genice2.cell import cellvectors
import genice2.lattices
import numpy as np
desc = {
    "ref": {
        "15_2_201714": "Engel 2018",
        "engel28": "Engel 2018"
    },
    "usage": "No options available.",
    "brief": "Hypothetical zeolitic ice"
}


class Lattice(genice2.lattices.Lattice):
    def __init__(self):
        self.cell = np.array([
            [1.192347, -1.336676, 3.388231],
            [0.882875, 6.031658, -0.189029],
            [-6.038371, -0.064357, 0.034045],
        ])
        self.waters = np.array([
            [0.202234, 0.470861, 0.178149],
            [0.344579, -0.054575, -0.330526],
            [-0.151382, -0.238174, -0.144889],
            [-0.293975, 0.287031, 0.363431],
            [0.267916, 0.352267, -0.232633],
            [-0.217342, -0.11977, 0.26581],
        ])
        self.coord = 'relative'
