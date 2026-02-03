from genice2.cell import cellvectors
import genice2.lattices
import numpy as np
desc = {
    "ref": {
        "12_2_32449": "Engel 2018",
        "engel09": "Engel 2018"
    },
    "usage": "No options available.",
    "brief": "Hypothetical zeolitic ice",
    "test": ({"args": "",
             "options": "-r 2 2 2"},)
}


class Lattice(genice2.lattices.Lattice):
    def __init__(self):
        self.cell = np.array([
            [4.906093, -2.166984, -1.768499],
            [-3.225085, 3.887306, -2.634777],
            [3.832468, 1.696101, 2.81583],
        ])
        self.waters = np.array([
            [0.156218, 0.39529, -0.383508],
            [-0.178719, -0.439869, 0.267773],
            [-0.445277, 0.020649, -0.324495],
            [0.431107, -0.051465, 0.202261],
            [0.226932, 0.248658, 0.010101],
            [-0.242547, -0.282667, -0.134704],
        ])
        self.coord = 'relative'
