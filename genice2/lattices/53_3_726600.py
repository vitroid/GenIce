desc={
    "ref": {
        "53_3_726600": "Engel, E.A., Anelli, A., Ceriotti, M. et al. Mapping uncharted territory in ice from zeolite networks to ice structures. Nat Commun 9, 2173 (2018). https://doi.org/10.1038/s41467-018-04618-6",
        "engel07": "Engel, E.A., Anelli, A., Ceriotti, M. et al. Mapping uncharted territory in ice from zeolite networks to ice structures. Nat Commun 9, 2173 (2018). https://doi.org/10.1038/s41467-018-04618-6"
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
            [13.982562, -4.048242, -0.028155],
            [-1.56494, 4.777418, 0.030278],
            [-0.00101, -0.03182, 4.704863],
        ])
        self.waters = np.array([
            [0.110619, -0.473564, 0.245892],
            [-0.029156, -0.239162, 0.25618],
            [0.388366, 0.480976, -0.251603],
            [-0.113279, 0.474237, -0.256748],
            [0.026141, 0.238496, -0.267066],
            [-0.388331, -0.469759, 0.244282],
            [0.472345, -0.230675, 0.235771],
            [-0.472598, 0.243982, -0.243614],
            [0.221838, 0.253899, 0.247924],
            [-0.278026, 0.253723, 0.245686],
            [0.276909, -0.247843, -0.252345],
            [-0.223715, -0.250359, -0.256199],
            [0.132511, -0.044199, -0.250777],
            [-0.134232, 0.045676, 0.243634],
            [0.3658, 0.050768, 0.249701],
            [-0.367525, -0.041869, -0.25449],
        ])
        self.coord = 'relative'
