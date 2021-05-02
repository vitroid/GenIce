from genice2.cell import cellvectors
import genice2.lattices
import numpy as np
desc = {
    "ref": {
        "169_2_7915": "Engel 2018",
        "engel06": "Engel 2018"
    },
    "usage": "No options available.",
    "brief": "Hypothetical zeolitic ice"
}


class Lattice(genice2.lattices.Lattice):
    def __init__(self):
        self.cell = np.array([
            [5.113806, 8.999167, 0.019509],
            [-5.25281, 9.003256, 0.009643],
            [-0.011739, -0.014666, 4.538722],
        ])
        self.waters = np.array([
            [-0.49089, -0.225094, 0.02519],
            [0.460436, 0.31797, -0.474697],
            [0.238915, 0.516278, -0.141744],
            [-0.27126, -0.425484, 0.36021],
            [-0.286816, -0.209004, -0.303936],
            [0.260365, 0.29656, 0.191742],
            [0.466271, 0.110763, -0.161188],
            [0.510635, -0.021299, 0.338924],
            [-0.47284, -0.430725, -0.32323],
            [0.439607, -0.475458, 0.17594],
            [0.051609, 0.495351, -0.491172],
            [-0.080399, -0.414873, 0.008093],
        ])
        self.coord = 'relative'
