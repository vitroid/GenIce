from genice2.cell import cellvectors
import genice2.lattices
import numpy as np
desc = {
    "ref": {
        "11_2_15848": "Engel 2018",
        "engel18": "Engel 2018"
    },
    "usage": "No options available.",
    "brief": "Hypothetical zeolitic ice",
    "test": ({"args": "",
             "options": "-r 2 2 2"},)
}


class Lattice(genice2.lattices.Lattice):
    def __init__(self):
        self.cell = np.array([
            [-2.319854, 0.93938, 5.267056],
            [0.83573, 5.812061, 1.143502],
            [-5.809818, -0.792117, 0.710229],
        ])
        self.waters = np.array([
            [-0.290188, 0.044868, 0.096462],
            [0.374418, -0.089897, 0.008098],
            [-0.327585, 0.507982, 0.244356],
            [0.439326, -0.557115, -0.12763],
            [0.229667, -0.067253, -0.337094],
            [-0.138403, 0.016385, 0.436419],
            [0.084123, 0.409276, -0.266745],
            [0.021352, -0.447365, 0.387203],
        ])
        self.coord = 'relative'
