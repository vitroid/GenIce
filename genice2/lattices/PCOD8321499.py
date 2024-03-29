from genice2.cell import cellvectors
import genice2.lattices
import numpy as np
desc = {
    "ref": {
        "PCOD8321499": "Engel 2018",
        "engel26": "Engel 2018"
    },
    "usage": "No options available.",
    "brief": "Hypothetical zeolitic ice"
}


class Lattice(genice2.lattices.Lattice):
    def __init__(self):
        self.cell = cellvectors(
            a=20.09717,
            b=20.04341,
            c=20.95821,
            A=89.5522,
            B=89.6358,
            C=60.9317
        )
        self.waters = np.array([
            [0.488392, -0.145329, 0.404387],
            [0.218323, -0.218308, 0.148320],
            [0.496688, -0.356575, 0.055932],
            [-0.171991, -0.015769, 0.325318],
            [0.209368, -0.011642, 0.308165],
            [-0.117154, 0.453008, 0.392823],
            [0.030348, 0.172379, 0.488752],
            [0.073805, 0.456589, 0.233702],
            [-0.124930, -0.358860, 0.074103],
            [0.492306, 0.453049, 0.392440],
            [0.216450, -0.409123, 0.314221],
            [0.679844, -0.148704, 0.230050],
            [0.487796, 0.047896, 0.069571],
            [0.632770, 0.179649, 0.481680],
            [0.214257, 0.183722, 0.151263],
            [0.687112, 0.453543, 0.229988],
            [-0.171827, 0.171851, 0.148306],
            [0.027318, -0.222809, 0.478221],
            [0.488737, -0.144275, -0.092670],
            [0.220967, -0.214751, 0.648680],
            [0.499913, -0.359152, 0.555487],
            [-0.172026, -0.015895, -0.178239],
            [0.211017, -0.008350, -0.191810],
            [-0.116996, 0.452866, -0.104109],
            [0.027599, 0.173436, -0.011312],
            [0.070334, 0.460757, 0.733576],
            [-0.127029, -0.354545, 0.573695],
            [0.492425, 0.452916, -0.104289],
            [0.218653, -0.409213, -0.189391],
            [0.677734, -0.149389, 0.730651],
            [0.484568, 0.046995, 0.569816],
            [0.636751, 0.177648, -0.018353],
            [0.211143, 0.187636, 0.651696],
            [0.686981, 0.453682, 0.733547],
            [-0.174047, 0.174696, 0.648026],
            [0.026550, -0.224306, -0.025459],
        ])
        self.coord = 'relative'
