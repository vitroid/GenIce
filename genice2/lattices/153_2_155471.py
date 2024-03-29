from genice2.cell import cellvectors
import genice2.lattices
import numpy as np
desc = {
    "ref": {
        "153_2_155471": "Engel 2018"
    },
    "usage": "No options available.",
    "brief": "Hypothetical zeolitic ice"
}


class Lattice(genice2.lattices.Lattice):
    def __init__(self):
        self.cell = cellvectors(
            a=17.71735,
            b=26.32421,
            c=18.61443,
            A=90.0,
            B=90.0,
            C=119.0366
        )
        self.waters = np.array([
            [0.008632, 0.166908, -0.083343],
            [0.254322, 0.331190, 0.253287],
            [0.250090, 0.166942, -0.416410],
            [0.499950, 0.166587, 0.416696],
            [0.255987, 0.003148, -0.246932],
            [0.258278, 0.166478, 0.083561],
            [0.508639, 0.166913, -0.083345],
            [-0.245684, 0.331187, 0.253288],
            [-0.249909, 0.166940, -0.416409],
            [-0.000052, 0.166585, 0.416697],
            [-0.244023, 0.003140, -0.246931],
            [-0.241722, 0.166478, 0.083562],
            [0.008632, 0.500241, -0.083348],
            [0.254326, 0.664526, 0.253286],
            [0.250091, 0.500276, -0.416409],
            [0.499948, 0.499915, 0.416692],
            [0.255986, 0.336480, -0.246933],
            [0.258275, 0.499806, 0.083560],
            [0.508629, 0.500241, -0.083347],
            [-0.245674, 0.664526, 0.253286],
            [-0.249908, 0.500276, -0.416410],
            [-0.000050, 0.499915, 0.416690],
            [-0.244017, 0.336479, -0.246939],
            [-0.241721, 0.499811, 0.083562],
            [0.008630, -0.166425, -0.083346],
            [0.254335, -0.002136, 0.253290],
            [0.250094, -0.166392, -0.416408],
            [0.499946, 0.833251, 0.416693],
            [0.255995, 0.669818, -0.246936],
            [0.258276, 0.833142, 0.083562],
            [0.508622, -0.166429, -0.083349],
            [-0.245680, -0.002144, 0.253290],
            [-0.249905, -0.166393, -0.416409],
            [-0.000048, 0.833255, 0.416701],
            [-0.244012, 0.669813, -0.246940],
            [-0.241723, 0.833142, 0.083560],
        ])
        self.coord = 'relative'
