desc={
    "ref": {
        "151_2_4949650": "Engel 2018",
        "engel23": "Engel 2018"
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
            a=22.74119,
            b=22.8581,
            c=23.48407,
            A=90.0,
            B=90.0,
            C=120.50374
        )
        self.waters = np.array([
            [0.453352, 0.120104, 0.346233],
            [0.160791, 0.043626, 0.176215],
            [0.072603, 0.286283, 0.060563],
            [0.213376, 0.286725, 0.227097],
            [0.207435, 0.426948, 0.392969],
            [0.146565, 0.106730, 0.460521],
            [0.455758, 0.353736, 0.296948],
            [0.380000, 0.339904, 0.012966],
            [0.393486, 0.043612, 0.130151],
            [-0.046649, 0.120104, 0.346234],
            [-0.339206, 0.043627, 0.176215],
            [0.572603, 0.286282, 0.060562],
            [-0.286624, 0.286724, 0.227098],
            [-0.292565, 0.426948, 0.392971],
            [-0.353435, 0.106729, 0.460521],
            [-0.044244, 0.353733, 0.296949],
            [-0.120000, 0.339903, 0.012966],
            [-0.106512, 0.043612, 0.130151],
            [0.453351, 0.620104, 0.346234],
            [0.160790, 0.543625, 0.176215],
            [0.072603, -0.213718, 0.060562],
            [0.213375, -0.213275, 0.227097],
            [0.207436, -0.073053, 0.392970],
            [0.146564, 0.606729, 0.460521],
            [0.455757, -0.146263, 0.296949],
            [0.380000, -0.160096, 0.012966],
            [0.393485, 0.543610, 0.130151],
            [-0.046648, 0.620104, 0.346234],
            [-0.339209, 0.543625, 0.176215],
            [0.572605, -0.213717, 0.060562],
            [-0.286624, -0.213275, 0.227096],
            [-0.292564, -0.073053, 0.392970],
            [-0.353435, 0.606729, 0.460521],
            [-0.044244, -0.146269, 0.296947],
            [-0.120000, -0.160096, 0.012966],
            [-0.106509, 0.543615, 0.130151],
            [0.453352, 0.120105, -0.153766],
            [0.160791, 0.043626, 0.676216],
            [0.072603, 0.286283, 0.560563],
            [0.213377, 0.286726, 0.727097],
            [0.207436, 0.426947, -0.107030],
            [0.146565, 0.106730, -0.039478],
            [0.455758, 0.353737, -0.203052],
            [0.380001, 0.339904, 0.512965],
            [0.393486, 0.043611, 0.630151],
            [-0.046649, 0.120103, -0.153767],
            [-0.339206, 0.043628, 0.676216],
            [0.572603, 0.286282, 0.560563],
            [-0.286623, 0.286726, 0.727097],
            [-0.292564, 0.426949, -0.107030],
            [-0.353435, 0.106730, -0.039478],
            [-0.044245, 0.353732, -0.203052],
            [-0.120000, 0.339903, 0.512965],
            [-0.106513, 0.043611, 0.630151],
            [0.453352, 0.620104, -0.153766],
            [0.160791, 0.543625, 0.676216],
            [0.072603, -0.213718, 0.560563],
            [0.213376, -0.213275, 0.727097],
            [0.207436, -0.073052, -0.107030],
            [0.146565, 0.606729, -0.039477],
            [0.455757, -0.146264, -0.203052],
            [0.380001, -0.160095, 0.512965],
            [0.393485, 0.543610, 0.630151],
            [-0.046648, 0.620104, -0.153766],
            [-0.339206, 0.543630, 0.676216],
            [0.572605, -0.213718, 0.560563],
            [-0.286624, -0.213274, 0.727097],
            [-0.292564, -0.073053, -0.107030],
            [-0.353435, 0.606729, -0.039478],
            [-0.044244, -0.146268, -0.203053],
            [-0.119999, -0.160097, 0.512965],
            [-0.106512, 0.543610, 0.630151],
        ])
        self.coord = 'relative'
