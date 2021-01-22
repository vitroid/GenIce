desc={
    "ref": {
        "PCOD8204698": "Engel 2018"
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
            a=17.2872,
            b=17.63854,
            c=44.42925,
            A=90.0,
            B=90.9037,
            C=90.0
        )
        self.waters = np.array([
            [0.004775, 0.245751, 0.082551],
            [0.124705, 0.400317, -0.166606],
            [0.375133, 0.433717, 0.329824],
            [0.124727, 0.099673, 0.333383],
            [0.495616, 0.229210, -0.085073],
            [0.375145, 0.066314, -0.170168],
            [0.380561, 0.088158, 0.168421],
            [0.380601, 0.411784, -0.331636],
            [0.250179, 0.241203, -0.250594],
            [0.250178, 0.258769, 0.249423],
            [0.133085, 0.424523, 0.163199],
            [0.133059, 0.075497, -0.336816],
            [0.175083, 0.116566, 0.000026],
            [0.341172, 0.133654, 0.496518],
            [0.337104, 0.363043, -0.002031],
            [0.175565, 0.382859, 0.500003],
            [0.004799, 0.253931, -0.417673],
            [0.496004, 0.271080, 0.414731],
            [0.504775, 0.245750, 0.082551],
            [0.624702, 0.400317, -0.166606],
            [-0.124867, 0.433718, 0.329824],
            [0.624723, 0.099673, 0.333383],
            [-0.004384, 0.229209, -0.085073],
            [-0.124856, 0.066314, -0.170168],
            [-0.119440, 0.088158, 0.168420],
            [-0.119397, 0.411784, -0.331636],
            [-0.249821, 0.241204, -0.250594],
            [-0.249821, 0.258770, 0.249423],
            [0.633085, 0.424523, 0.163200],
            [0.633057, 0.075497, -0.336816],
            [0.675084, 0.116566, 0.000026],
            [-0.158829, 0.133656, 0.496518],
            [-0.162896, 0.363044, -0.002031],
            [0.675567, 0.382860, 0.500003],
            [0.504800, 0.253931, -0.417673],
            [-0.003996, 0.271080, 0.414731],
            [0.004774, 0.745748, 0.082551],
            [0.124704, -0.099683, -0.166606],
            [0.375133, -0.066283, 0.329824],
            [0.124726, 0.599675, 0.333383],
            [0.495615, 0.729210, -0.085073],
            [0.375144, 0.566313, -0.170168],
            [0.380560, 0.588155, 0.168420],
            [0.380603, -0.088216, -0.331636],
            [0.250179, 0.741201, -0.250594],
            [0.250179, -0.241230, 0.249423],
            [0.133085, -0.075477, 0.163200],
            [0.133059, 0.575495, -0.336816],
            [0.175082, 0.616565, 0.000026],
            [0.341172, 0.633652, 0.496518],
            [0.337103, -0.136956, -0.002031],
            [0.175566, -0.117141, 0.500003],
            [0.004799, -0.246068, -0.417673],
            [0.496004, -0.228920, 0.414731],
            [0.504775, 0.745753, 0.082551],
            [0.624702, -0.099683, -0.166606],
            [-0.124866, -0.066283, 0.329824],
            [0.624729, 0.599675, 0.333383],
            [-0.004384, 0.729210, -0.085073],
            [-0.124857, 0.566312, -0.170169],
            [-0.119440, 0.588155, 0.168421],
            [-0.119398, -0.088215, -0.331636],
            [-0.249821, 0.741201, -0.250594],
            [-0.249822, -0.241231, 0.249423],
            [0.633085, -0.075478, 0.163200],
            [0.633057, 0.575495, -0.336816],
            [0.675084, 0.616565, 0.000026],
            [-0.158828, 0.633652, 0.496518],
            [-0.162896, -0.136956, -0.002031],
            [0.675567, -0.117140, 0.500003],
            [0.504799, -0.246068, -0.417673],
            [-0.003996, -0.228920, 0.414731],
        ])
        self.coord = 'relative'
