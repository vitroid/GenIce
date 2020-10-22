desc={
    "ref": {
        "PCOD8047931": "Engel, E.A., Anelli, A., Ceriotti, M. et al. Mapping uncharted territory in ice from zeolite networks to ice structures. Nat Commun 9, 2173 (2018). https://doi.org/10.1038/s41467-018-04618-6",
        "engel27": "Engel, E.A., Anelli, A., Ceriotti, M. et al. Mapping uncharted territory in ice from zeolite networks to ice structures. Nat Commun 9, 2173 (2018). https://doi.org/10.1038/s41467-018-04618-6"
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
            [3.066949, 0.001744, -5.534685],
            [0.004297, 6.987767, 0.007504],
            [7.2388, -0.001, 5.255684],
        ])
        self.waters = np.array([
            [0.120267, -0.003589, 0.413046],
            [-0.09487, 0.497214, -0.3731],
            [-0.445567, -0.011575, 0.12533],
            [0.470605, -0.510982, -0.087173],
            [-0.176505, 0.173521, -0.484494],
            [0.201653, -0.328522, -0.478157],
            [-0.278708, -0.005907, -0.259514],
            [0.302557, 0.493308, 0.296477],
            [0.074373, -0.429759, -0.104851],
            [-0.048737, 0.069011, 0.144133],
            [-0.436445, -0.334207, -0.322604],
            [0.461107, 0.165695, 0.360369],
            [-0.051607, 0.435024, 0.131946],
            [0.075423, -0.063589, -0.094076],
            [-0.38308, -0.372899, 0.176824],
            [0.407938, 0.127952, -0.13876],
        ])
        self.coord = 'relative'