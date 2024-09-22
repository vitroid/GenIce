from genice2.cell import cellvectors
import genice2.lattices
import numpy as np
desc = {
    "ref": {
        "PCOD8047078": "Engel 2018",
        "engel13": "Engel 2018"
    },
    "usage": "No options available.",
    "brief": "Hypothetical zeolitic ice"
}


class Lattice(genice2.lattices.Lattice):
    def __init__(self):
        self.cell = np.array([
            [1.871528, 0.031133, -3.760368],
            [-0.103113, 9.627472, 0.029368],
            [7.102959, -0.011856, 3.350493],
        ])
        self.waters = np.array([
            [-0.447084, -0.483599, 0.138541],
            [0.442916, 0.029262, -0.157831],
            [0.49811, 0.005808, 0.159092],
            [-0.509321, -0.481791, -0.178883],
            [-0.467688, -0.240919, 0.260863],
            [0.461353, 0.269874, -0.2807],
            [0.024071, -0.228034, 0.44058],
            [-0.03107, 0.272968, 0.538244],
            [-0.005907, -0.357623, -0.27061],
            [-0.00587, 0.128645, 0.255657],
            [-0.051355, -0.084965, -0.275396],
            [0.046943, 0.400465, 0.247554],
        ])
        self.coord = 'relative'
