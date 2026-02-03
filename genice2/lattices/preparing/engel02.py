from genice2.cell import cellvectors
import genice2.lattices
import numpy as np
desc = {
    "ref": {
        "12_2_29187": "Engel 2018",
        "engel02": "Engel 2018"
    },
    "usage": "No options available.",
    "brief": "Hypothetical zeolitic ice",
    "test": ({"args": "",
             "options": "-r 2 2 2 --depol=none"},)
}


class Lattice(genice2.lattices.Lattice):
    def __init__(self):
        self.cell = np.array([
            [6.151316, -0.027655, -0.035754],
            [-2.042472, 5.845221, 0.026555],
            [-1.978529, -2.930206, 5.013556],
        ])
        self.waters = np.array([
            [-0.011996, 0.225767, -0.246008],
            [-0.02371, -0.272108, 0.24034],
            [-0.322387, -0.068408, 0.250959],
            [0.2864, 0.023605, -0.255622],
            [0.182996, -0.274515, -0.058655],
            [-0.21879, 0.229406, 0.054395],
        ])
        self.coord = 'relative'
