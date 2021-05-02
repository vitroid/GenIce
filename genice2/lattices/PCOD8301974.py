from genice2.cell import cellvectors
import genice2.lattices
import numpy as np
desc = {
    "ref": {
        "PCOD8301974": "Engel 2018",
        "engel20": "Engel 2018"
    },
    "usage": "No options available.",
    "brief": "Hypothetical zeolitic ice"
}


class Lattice(genice2.lattices.Lattice):
    def __init__(self):
        self.cell = cellvectors(
            a=30.61955,
            b=30.70538,
            c=24.56766,
            A=90.0,
            B=90.0,
            C=90.0
        )
        self.waters = np.array([
            [0.292711, 0.327630, 0.328359],
            [0.077596, 0.447710, 0.458224],
            [0.045100, 0.084941, 0.081428],
            [0.442701, 0.412882, 0.333948],
            [0.162623, 0.297621, 0.453303],
            [0.042475, 0.332274, 0.190743],
            [0.297235, 0.082721, 0.440437],
            [0.196674, 0.167589, 0.078141],
            [0.325163, 0.449807, 0.062846],
            [0.192852, 0.413257, 0.190851],
            [0.447737, 0.167676, 0.440306],
            [0.162895, 0.047346, 0.315923],
            [0.079667, 0.199670, 0.312666],
            [0.412113, 0.302224, 0.065617],
            [0.411995, 0.046952, 0.208659],
            [0.328071, 0.197764, 0.202870],
            [-0.207289, 0.327630, 0.328359],
            [0.577595, 0.447710, 0.458224],
            [0.545099, 0.084941, 0.081428],
            [-0.057298, 0.412882, 0.333948],
            [0.662622, 0.297621, 0.453303],
            [0.542474, 0.332274, 0.190743],
            [-0.202765, 0.082721, 0.440437],
            [0.696673, 0.167589, 0.078141],
            [-0.174836, 0.449807, 0.062847],
            [0.692851, 0.413257, 0.190851],
            [-0.052263, 0.167676, 0.440306],
            [0.662897, 0.047346, 0.315923],
            [0.579666, 0.199670, 0.312666],
            [-0.087888, 0.302224, 0.065617],
            [-0.088004, 0.046952, 0.208659],
            [-0.171928, 0.197764, 0.202870],
            [0.292708, -0.172367, 0.328361],
            [0.077604, -0.052262, 0.458192],
            [0.045104, 0.584943, 0.081446],
            [0.442698, -0.087115, 0.333953],
            [0.162616, -0.202377, 0.453303],
            [0.042459, -0.167732, 0.190750],
            [0.297224, 0.582722, 0.440449],
            [0.196583, 0.667561, 0.078139],
            [0.325171, -0.050182, 0.062852],
            [0.192873, -0.086779, 0.190833],
            [0.447747, 0.667681, 0.440311],
            [0.162891, 0.547354, 0.315940],
            [0.079667, 0.699669, 0.312662],
            [0.412099, -0.197786, 0.065621],
            [0.411972, 0.546920, 0.208629],
            [0.328062, 0.697777, 0.202850],
            [-0.207292, -0.172367, 0.328361],
            [0.577605, -0.052262, 0.458192],
            [0.545103, 0.584943, 0.081446],
            [-0.057301, -0.087116, 0.333953],
            [0.662616, -0.202377, 0.453303],
            [0.542457, -0.167732, 0.190749],
            [-0.202776, 0.582722, 0.440453],
            [0.696584, 0.667561, 0.078138],
            [-0.174829, -0.050182, 0.062852],
            [0.692871, -0.086779, 0.190833],
            [-0.052252, 0.667681, 0.440311],
            [0.662890, 0.547354, 0.315939],
            [0.579669, 0.699669, 0.312662],
            [-0.087899, -0.197786, 0.065621],
            [-0.088028, 0.546920, 0.208629],
            [-0.171939, 0.697777, 0.202850],
            [0.293072, 0.328206, -0.171776],
            [0.077828, 0.452315, -0.041904],
            [0.047197, 0.082519, 0.578846],
            [0.444755, 0.415155, -0.168605],
            [0.164685, 0.299925, -0.043872],
            [0.042770, 0.332385, 0.690326],
            [0.297194, 0.082675, -0.059486],
            [0.196680, 0.167416, 0.578049],
            [0.325052, 0.449677, 0.562626],
            [0.192910, 0.416904, 0.690611],
            [0.447789, 0.167818, -0.059972],
            [0.163105, 0.048110, -0.183851],
            [0.079925, 0.200025, -0.186880],
            [0.410091, 0.300100, 0.568463],
            [0.412472, 0.047113, 0.703449],
            [0.327693, 0.197448, 0.703421],
            [-0.206927, 0.328203, -0.171775],
            [0.577827, 0.452315, -0.041905],
            [0.547196, 0.082519, 0.578846],
            [-0.055245, 0.415155, -0.168605],
            [0.664683, 0.299925, -0.043872],
            [0.542771, 0.332381, 0.690330],
            [-0.202806, 0.082675, -0.059486],
            [0.696679, 0.167416, 0.578049],
            [-0.174947, 0.449677, 0.562626],
            [0.692914, 0.416901, 0.690615],
            [-0.052212, 0.167818, -0.059972],
            [0.663106, 0.048110, -0.183851],
            [0.579924, 0.200025, -0.186879],
            [-0.089911, 0.300100, 0.568463],
            [-0.087527, 0.047113, 0.703449],
            [-0.172308, 0.197448, 0.703425],
            [0.292490, -0.172166, -0.171740],
            [0.077337, -0.048252, -0.041284],
            [0.047210, 0.582553, 0.578814],
            [0.442760, -0.086310, -0.166615],
            [0.162464, -0.201871, -0.046667],
            [0.042518, -0.167579, 0.690595],
            [0.297310, 0.582742, -0.059432],
            [0.196706, 0.667424, 0.578036],
            [0.325049, -0.050316, 0.562630],
            [0.192903, -0.083106, 0.690668],
            [0.447740, 0.667824, -0.059951],
            [0.165081, 0.549923, -0.187034],
            [0.079722, 0.699855, -0.187341],
            [0.410078, -0.199884, 0.568483],
            [0.412400, 0.547341, 0.703201],
            [0.327673, 0.697432, 0.703384],
            [-0.207511, -0.172165, -0.171741],
            [0.577337, -0.048252, -0.041284],
            [0.547209, 0.582553, 0.578814],
            [-0.057240, -0.086310, -0.166615],
            [0.662466, -0.201871, -0.046667],
            [0.542519, -0.167579, 0.690595],
            [-0.202690, 0.582742, -0.059432],
            [0.696705, 0.667424, 0.578036],
            [-0.174951, -0.050316, 0.562630],
            [0.692904, -0.083106, 0.690668],
            [-0.052261, 0.667824, -0.059951],
            [0.665082, 0.549920, -0.187033],
            [0.579721, 0.699855, -0.187341],
            [-0.089921, -0.199884, 0.568483],
            [-0.087601, 0.547341, 0.703201],
            [-0.172328, 0.697432, 0.703384],
        ])
        self.coord = 'relative'
