from genice2.cell import cellvectors
import genice2.lattices
import numpy as np
desc = {
    "ref": {
        "PCOD8045578": "Engel 2018",
        "engel21": "Engel 2018"
    },
    "usage": "No options available.",
    "brief": "Hypothetical zeolitic ice",
    "test": ({"args": "",
             "options": "-r 2 2 2"},)
}


class Lattice(genice2.lattices.Lattice):
    def __init__(self):
        self.cell = np.array([
            [7.734411, -0.004397, -1.982279],
            [-0.00197, 4.250365, 0.000625],
            [6.795719, -0.000229, 3.376358],
        ])
        self.waters = np.array([
            [0.534822, 0.377526, -0.311948],
            [-0.565464, 0.37787, 0.309548],
            [-0.065053, -0.11937, 0.308921],
            [0.036048, -0.119884, -0.313818],
            [-0.322001, -0.119833, 0.221843],
            [0.292833, -0.121186, -0.226351],
            [0.175819, 0.379229, 0.224696],
            [-0.207129, 0.378088, -0.226363],
        ])
        self.coord = 'relative'
