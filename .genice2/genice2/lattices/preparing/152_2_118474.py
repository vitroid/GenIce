from genice2.cell import cellvectors
import genice2.lattices
import numpy as np
desc = {
    "ref": {
        "152_2_118474": "Engel 2018",
        "engel16": "Engel 2018"
    },
    "usage": "No options available.",
    "brief": "Hypothetical zeolitic ice"
}


class Lattice(genice2.lattices.Lattice):
    def __init__(self):
        self.cell = np.array([
            [2.749962, 5.015182, 0.142688],
            [-2.96902, 4.885153, 0.280579],
            [0.231039, -0.302749, 5.926585],
        ])
        self.waters = np.array([
            [-0.39041, -0.13446, 0.2867],
            [0.327705, 0.513437, 0.083761],
            [0.112064, 0.383286, 0.413665],
            [-0.522271, -0.350907, -0.382067],
            [-0.175229, -0.484807, -0.051938],
            [0.459917, 0.165008, -0.254772],
            [0.20679, 0.01487, 0.014993],
            [-0.248653, 0.237619, 0.345489],
            [-0.027272, -0.20832, -0.30977],
        ])
        self.coord = 'relative'
