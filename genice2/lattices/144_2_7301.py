from genice2.cell import cellvectors
import genice2.lattices
import numpy as np
desc = {
    "ref": {
        "144_2_7301": "Engel 2018"
    },
    "usage": "No options available.",
    "brief": "Hypothetical zeolitic ice"
}


class Lattice(genice2.lattices.Lattice):
    def __init__(self):
        self.cell = cellvectors(
            a=23.78803,
            b=23.78803,
            c=23.01145,
            A=90.0,
            B=90.0,
            C=120.0
        )
        self.waters = np.array([
            [0.381880, 0.115750, 0.089262],
            [0.383686, 0.265280, 0.255759],
            [0.234555, 0.118136, 0.422564],
            [0.118700, 0.381677, 0.340143],
            [0.262976, 0.381255, 0.173517],
            [0.118289, 0.236973, 0.006808],
            [-0.118120, 0.115751, 0.089262],
            [-0.116314, 0.265281, 0.255759],
            [-0.265445, 0.118136, 0.422564],
            [0.618700, 0.381676, 0.340143],
            [-0.237025, 0.381255, 0.173516],
            [0.618288, 0.236973, 0.006808],
            [0.381880, 0.615751, 0.089262],
            [0.383686, -0.234720, 0.255758],
            [0.234556, 0.618134, 0.422564],
            [0.118700, -0.118323, 0.340143],
            [0.262977, -0.118745, 0.173516],
            [0.118289, 0.736973, 0.006808],
            [-0.118121, 0.615751, 0.089262],
            [-0.116315, -0.234721, 0.255758],
            [-0.265445, 0.618134, 0.422564],
            [0.618700, -0.118323, 0.340143],
            [-0.237024, -0.118745, 0.173516],
            [0.618289, 0.736973, 0.006808],
            [0.381880, 0.115751, 0.589263],
            [0.383686, 0.265281, -0.244241],
            [0.234555, 0.118136, -0.077436],
            [0.118700, 0.381677, -0.159857],
            [0.262976, 0.381255, 0.673517],
            [0.118289, 0.236973, 0.506809],
            [-0.118120, 0.115750, 0.589263],
            [-0.116314, 0.265280, -0.244241],
            [-0.265445, 0.118136, -0.077436],
            [0.618700, 0.381676, -0.159857],
            [-0.237025, 0.381255, 0.673517],
            [0.618288, 0.236973, 0.506809],
            [0.381881, 0.615751, 0.589263],
            [0.383685, -0.234722, -0.244242],
            [0.234555, 0.618134, -0.077437],
            [0.118700, -0.118324, -0.159857],
            [0.262977, -0.118745, 0.673517],
            [0.118289, 0.736973, 0.506809],
            [-0.118121, 0.615751, 0.589263],
            [-0.116316, -0.234722, -0.244243],
            [-0.265445, 0.618134, -0.077436],
            [0.618700, -0.118324, -0.159857],
            [-0.237024, -0.118745, 0.673517],
            [0.618289, 0.736973, 0.506809],
        ])
        self.coord = 'relative'
