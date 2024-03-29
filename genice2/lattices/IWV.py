from genice2.cell import cellvectors
import genice2.lattices
import numpy as np
desc = {
    "ref": {
        "IWV": "Engel 2018"
    },
    "usage": "No options available.",
    "brief": "Hypothetical zeolitic ice"
}


class Lattice(genice2.lattices.Lattice):
    def __init__(self):
        self.cell = cellvectors(
            a=23.472,
            b=24.65765,
            c=26.95031,
            A=77.3042,
            B=63.4086,
            C=61.4293
        )
        self.waters = np.array([
            [0.392703, -0.199826, -0.083268],
            [0.367098, 0.187920, 0.194027],
            [0.219961, 0.381152, 0.100095],
            [0.019925, -0.000255, -0.105494],
            [0.620319, -0.131433, -0.389438],
            [0.558508, -0.192507, 0.194635],
            [0.986497, -0.000164, 0.100291],
            [0.499887, 0.112932, -0.394491],
            [0.419588, -0.000414, 0.095047],
            [0.593512, 0.000156, -0.111290],
            [0.613213, 0.112538, 0.383597],
            [0.939825, -0.198423, 0.193457],
            [0.233226, 0.000071, 0.038928],
            [0.793920, 0.386971, 0.099873],
            [0.779410, -0.005963, -0.049479],
            [0.213169, -0.387527, -0.111004],
            [0.446734, 0.193706, -0.205426],
            [0.613573, 0.192327, 0.078395],
            [0.633656, -0.193502, -0.194960],
            [0.826965, 0.387251, -0.111187],
            [0.732260, -0.130545, 0.389678],
            [0.039550, 0.387773, 0.033362],
            [0.406694, 0.381255, -0.110972],
            [0.506531, -0.125035, 0.388652],
            [0.967135, -0.387875, -0.050029],
            [0.273485, 0.118614, -0.394362],
            [0.419887, -0.393473, 0.044763],
            [0.393192, -0.125695, -0.393519],
            [0.254010, -0.193509, -0.194826],
            [0.579895, 0.387964, -0.044961],
            [0.600039, -0.381556, 0.104996],
            [0.200008, 0.193706, -0.094156],
            [0.746856, 0.193150, 0.189364],
            [0.380568, 0.118551, 0.388468],
            [0.807491, -0.199884, 0.082959],
            [0.779882, -0.381614, -0.105355],
            [0.174832, -0.388224, 0.099095],
            [0.065757, 0.188842, -0.200169],
        ])
        self.coord = 'relative'
