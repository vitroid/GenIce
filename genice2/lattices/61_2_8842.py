desc={
    "ref": {
        "61_2_8842": "Engel, E.A., Anelli, A., Ceriotti, M. et al. Mapping uncharted territory in ice from zeolite networks to ice structures. Nat Commun 9, 2173 (2018). https://doi.org/10.1038/s41467-018-04618-6",
        "engel11": "Engel, E.A., Anelli, A., Ceriotti, M. et al. Mapping uncharted territory in ice from zeolite networks to ice structures. Nat Commun 9, 2173 (2018). https://doi.org/10.1038/s41467-018-04618-6"
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
            [7.125918, -0.079424, 0.058654],
            [0.034373, 7.104955, -0.037821],
            [0.008189, -0.216, 7.058855],
        ])
        self.waters = np.array([
            [0.109002, 0.106928, 0.113962],
            [0.422775, -0.084517, -0.416623],
            [-0.097484, -0.116811, -0.090457],
            [0.40433, -0.398191, 0.083551],
            [0.104564, 0.377627, -0.426758],
            [-0.364147, 0.134528, 0.388068],
            [-0.383985, 0.389896, -0.129171],
            [-0.095589, -0.40074, 0.380522],
            [0.406747, 0.408518, 0.380201],
            [0.126478, -0.395123, -0.140986],
            [-0.384293, -0.385304, -0.402943],
            [0.125981, -0.117373, 0.385678],
            [0.390208, 0.112263, -0.109423],
            [-0.086728, 0.401981, 0.070029],
            [-0.090506, 0.076408, -0.388842],
            [-0.404667, -0.090822, 0.100596],
        ])
        self.coord = 'relative'
