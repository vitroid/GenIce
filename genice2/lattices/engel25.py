desc={
    "ref": {
        "2_2_342692": "Engel, E.A., Anelli, A., Ceriotti, M. et al. Mapping uncharted territory in ice from zeolite networks to ice structures. Nat Commun 9, 2173 (2018). https://doi.org/10.1038/s41467-018-04618-6",
        "engel25": "Engel, E.A., Anelli, A., Ceriotti, M. et al. Mapping uncharted territory in ice from zeolite networks to ice structures. Nat Commun 9, 2173 (2018). https://doi.org/10.1038/s41467-018-04618-6"
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
