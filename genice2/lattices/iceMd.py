#!/usr/bin/python
# coding: utf-8


from genice2.lattices.iceM import Lattice as LatticeM

desc = {
    "ref": {"Md": "Mochizuki 2024"},
    "usage": "No options available.",
    "brief": "A hydrogen-disordered counterpart of ice M.",
    "test": ({"options": "--depol=none"},),
}


class Lattice(LatticeM):
    def __init__(self):
        super().__init__()
        del self.fixed
