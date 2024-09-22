from genice2.cell import cellvectors
import genice2.lattices
import numpy as np
desc = {
    "ref": {
        "2_2_342692": "Engel 2018",
        "engel25": "Engel 2018"
    },
    "usage": "No options available.",
    "brief": "Hypothetical zeolitic ice"
}


class Lattice(genice2.lattices.Lattice):
    def __init__(self):
        self.cell = np.array([
            [0.368881, -2.911415, 0.057559],
            [1.833515, 0.322898, 4.525229],
            [-4.548334, -0.55816, 1.889219],
        ])
        self.waters = np.array([
            [-0.022226, 0.141384, 0.24748],
            [0.125021, -0.017838, -0.252682],
            [0.529194, -0.522451, -0.253412],
            [-0.426467, -0.354032, 0.246384],
        ])
        self.coord = 'relative'
