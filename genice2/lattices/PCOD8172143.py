from genice2.cell import cellvectors
import genice2.lattices
import numpy as np
desc = {
    "ref": {
        "PCOD8172143": "Engel 2018",
        "engel15": "Engel 2018"
    },
    "usage": "No options available.",
    "brief": "Hypothetical zeolitic ice"
}


class Lattice(genice2.lattices.Lattice):
    def __init__(self):
        self.cell = np.array([
            [11.008881, 0.009294, -0.003651],
            [-0.004573, 3.991893, 0.012244],
            [0.004025, -0.015182, 5.062874],
        ])
        self.waters = np.array([
            [-0.146857, -0.164981, -0.145346],
            [0.138435, 0.407606, -0.141405],
            [-0.362197, 0.325872, 0.354342],
            [0.352971, -0.095426, 0.359066],
            [-0.004425, 0.118257, 0.180275],
            [0.495285, -0.388743, -0.320714],
            [-0.364231, -0.170917, 0.039851],
            [0.355472, 0.401272, 0.044282],
            [-0.145002, 0.338687, -0.460097],
            [0.135439, -0.095384, -0.456185],
        ])
        self.coord = 'relative'
