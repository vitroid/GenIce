desc={
    "ref": {
        "2_2_623457": "Engel, E.A., Anelli, A., Ceriotti, M. et al. Mapping uncharted territory in ice from zeolite networks to ice structures. Nat Commun 9, 2173 (2018). https://doi.org/10.1038/s41467-018-04618-6"
    },
    "usage": "No options available.",
    "brief": "Hypothetical zeolitic ice"
}
import numpy as np
import genice2.lattices
from genice2.cell import cellvectors

class Lattice(genice2.lattices.Lattice):
    def __init__(self):
        self.cell = cellvectors(
            a=25.83251,
            b=29.04093,
            c=27.40942,
            A=90.0,
            B=90.0,
            C=90.0
        )
        self.waters = np.array([
            [0.166953, 0.083210, 0.011201],
            [0.163821, 0.418768, 0.263662],
            [0.000289, 0.333207, 0.011200],
            [-0.002839, 0.168759, 0.263656],
            [-0.003007, 0.336496, 0.324382],
            [0.001657, 0.164088, 0.074798],
            [0.163672, 0.086511, 0.324369],
            [0.168306, 0.414109, 0.074784],
            [0.500284, 0.083210, 0.011201],
            [0.497153, 0.418768, 0.263662],
            [0.333622, 0.333207, 0.011200],
            [0.330494, 0.168760, 0.263655],
            [0.330327, 0.336496, 0.324382],
            [0.334990, 0.164088, 0.074798],
            [0.497006, 0.086512, 0.324368],
            [0.501639, 0.414109, 0.074783],
            [-0.166381, 0.083210, 0.011201],
            [-0.169513, 0.418768, 0.263662],
            [0.666954, 0.333207, 0.011200],
            [0.663826, 0.168759, 0.263655],
            [0.663660, 0.336496, 0.324382],
            [0.668325, 0.164088, 0.074799],
            [-0.169661, 0.086511, 0.324368],
            [-0.165028, 0.414109, 0.074784],
            [0.166953, 0.583211, 0.011201],
            [0.163821, -0.081232, 0.263662],
            [0.000289, -0.166793, 0.011200],
            [-0.002839, 0.668760, 0.263655],
            [-0.003006, -0.163503, 0.324382],
            [0.001657, 0.664087, 0.074798],
            [0.163672, 0.586510, 0.324368],
            [0.168306, -0.085893, 0.074783],
            [0.500284, 0.583211, 0.011201],
            [0.497153, -0.081232, 0.263662],
            [0.333622, -0.166793, 0.011200],
            [0.330494, 0.668760, 0.263656],
            [0.330327, -0.163504, 0.324382],
            [0.334990, 0.664087, 0.074798],
            [0.497006, 0.586510, 0.324368],
            [0.501639, -0.085893, 0.074783],
            [-0.166381, 0.583211, 0.011201],
            [-0.169513, -0.081232, 0.263662],
            [0.666954, -0.166793, 0.011200],
            [0.663826, 0.668760, 0.263655],
            [0.663660, -0.163503, 0.324382],
            [0.668325, 0.664087, 0.074799],
            [-0.169661, 0.586510, 0.324368],
            [-0.165028, -0.085893, 0.074784],
            [0.161138, 0.083694, 0.511426],
            [0.172129, 0.416268, -0.239019],
            [-0.005527, 0.333700, 0.511441],
            [0.005476, 0.166273, -0.239024],
            [0.005465, 0.333667, -0.177656],
            [-0.005502, 0.166391, 0.572270],
            [0.172145, 0.083677, -0.177647],
            [0.161136, 0.416402, 0.572267],
            [0.494470, 0.083693, 0.511426],
            [0.505464, 0.416268, -0.239019],
            [0.327807, 0.333700, 0.511441],
            [0.338810, 0.166274, -0.239024],
            [0.338798, 0.333667, -0.177656],
            [0.327831, 0.166392, 0.572270],
            [0.505479, 0.083677, -0.177647],
            [0.494470, 0.416402, 0.572267],
            [-0.172195, 0.083693, 0.511426],
            [-0.161205, 0.416268, -0.239019],
            [0.661140, 0.333700, 0.511441],
            [0.672141, 0.166273, -0.239024],
            [0.672130, 0.333667, -0.177656],
            [0.661163, 0.166390, 0.572270],
            [-0.161189, 0.083677, -0.177647],
            [-0.172197, 0.416402, 0.572267],
            [0.161137, 0.583693, 0.511426],
            [0.172128, -0.083732, -0.239019],
            [-0.005527, -0.166300, 0.511441],
            [0.005477, 0.666273, -0.239024],
            [0.005465, -0.166333, -0.177656],
            [-0.005502, 0.666391, 0.572270],
            [0.172145, 0.583676, -0.177647],
            [0.161136, -0.083599, 0.572267],
            [0.494470, 0.583693, 0.511426],
            [0.505464, -0.083731, -0.239019],
            [0.327807, -0.166300, 0.511441],
            [0.338810, 0.666273, -0.239024],
            [0.338798, -0.166333, -0.177656],
            [0.327831, 0.666391, 0.572270],
            [0.505479, 0.583676, -0.177647],
            [0.494470, -0.083599, 0.572267],
            [-0.172196, 0.583693, 0.511426],
            [-0.161204, -0.083731, -0.239019],
            [0.661140, -0.166300, 0.511441],
            [0.672141, 0.666273, -0.239024],
            [0.672130, -0.166333, -0.177656],
            [0.661163, 0.666391, 0.572270],
            [-0.161189, 0.583676, -0.177647],
            [-0.172197, -0.083599, 0.572267],
        ])
        self.coord = 'relative'
