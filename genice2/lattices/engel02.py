desc={
    "ref": {
        "12_2_29187": "Engel, E.A., Anelli, A., Ceriotti, M. et al. Mapping uncharted territory in ice from zeolite networks to ice structures. Nat Commun 9, 2173 (2018). https://doi.org/10.1038/s41467-018-04618-6",
        "engel02": "Engel, E.A., Anelli, A., Ceriotti, M. et al. Mapping uncharted territory in ice from zeolite networks to ice structures. Nat Commun 9, 2173 (2018). https://doi.org/10.1038/s41467-018-04618-6"
    },
    "usage": "No options available.",
    "brief": "Hypothetical zeolitic ice"
}
import numpy as np
import genice2.lattices
from genice2.cell import cellvectors

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
