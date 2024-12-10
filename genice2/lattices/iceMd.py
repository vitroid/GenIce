#!/usr/bin/python
# coding: utf-8


import genice2

desc = {
    "ref": {"Md": "Mochizuki 2024"},
    "usage": "No options available.",
    "brief": "A hydrogen-disordered counterpart of ice M.",
    "test": ({"options": "--depol=none"},),
}


class Lattice(genice2.lattices.iceM.Lattice):
    def __init__(self):
        super().__init__()
        del self.fixed
