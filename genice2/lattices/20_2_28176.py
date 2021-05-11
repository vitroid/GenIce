from genice2.cell import cellvectors
import genice2.lattices
import numpy as np
desc = {
    "ref": {
        "20_2_28176": "Engel 2018",
        "engel32": "Engel 2018"
    },
    "usage": "No options available.",
    "brief": "Hypothetical zeolitic ice"
}


class Lattice(genice2.lattices.Lattice):
    def __init__(self):
        self.cell = np.array([
            [3.418243, 5.235012, -0.004297],
            [-3.428522, 5.250172, -0.022436],
            [0.000972, 0.362817, 3.96255],
        ])
        self.waters = np.array([
            [-7.2e-05, 0.358929, 0.130279],
            [-0.00188, -0.38112, -0.396162],
            [-0.399506, 0.37382, 0.345636],
            [0.368018, -0.393184, -0.153058],
            [-0.383528, -0.024985, 0.559425],
            [0.355719, -0.022976, 0.08705],
        ])
        self.coord = 'relative'
