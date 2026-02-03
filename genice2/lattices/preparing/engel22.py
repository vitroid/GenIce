from genice2.cell import cellvectors
import genice2.lattices
import numpy as np
desc = {
    "ref": {
        "58_2_511": "Engel 2018",
        "engel22": "Engel 2018"
    },
    "usage": "No options available.",
    "brief": "Hypothetical zeolitic ice"
}


class Lattice(genice2.lattices.Lattice):
    def __init__(self):
        self.cell = np.array([
            [12.596398, 0.051273, 0.037264],
            [-0.010725, 4.049569, -0.015741],
            [0.020623, 0.018489, 6.329635],
        ])
        self.waters = np.array([
            [0.183283, 0.247762, 0.20573],
            [-0.20373, -0.273234, -0.205554],
            [0.299427, -0.251977, 0.304239],
            [-0.199261, -0.279293, 0.195541],
            [0.296815, -0.239898, -0.293499],
            [-0.31835, 0.226834, -0.308766],
            [0.179605, 0.261641, -0.196397],
            [-0.317091, 0.222044, 0.289279],
            [0.047193, 0.237493, 0.504037],
            [-0.067604, -0.261368, 0.4973],
            [0.43262, -0.25698, 0.00197],
            [-0.451901, 0.242393, -0.007378],
        ])
        self.coord = 'relative'
