"""

Command line: /Users/matto/miniforge3/bin/genice2 xdtc2 -f reshape
Reshaping the unit cell.
  i:[1 0 0]
  j:[0 1 0]
  k:[0 0 1]
"""

import genice2.lattices
desc = {"ref": {"xdtc": 'Matsumoto 2021'},
        "usage": "No options available.",
        "brief": "A porous ice with cylindrical channels."
        }


class Lattice(genice2.lattices.Lattice):
    def __init__(self):
        self.bondlen = 0.2810520357047384
        self.coord = 'relative'
        from genice2.cell import cellvectors
        self.cell = cellvectors(a=2.04622653, b=3.54416832, c=0.70192855)
        self.density = 0.610666558384
        self.waters = """
    0.6023    0.1326    0.9343
    0.3977    0.1326    0.9343
    0.5000    0.2348    0.9343
    0.3977    0.1326    0.5657
    0.6023    0.1326    0.5657
    0.5000    0.2348    0.5657
    0.5000    0.1666    0.4343
    0.3207    0.1843    0.4343
    0.5631    0.0682    0.4343
    0.6161    0.2475    0.4343
    0.2816    0.2815    0.9343
    0.7815    0.2184    0.9343
    0.6793    0.1843    0.4343
    0.3839    0.2475    0.4343
    0.4368    0.0682    0.4343
    0.2185    0.2184    0.9343
    0.7184    0.2815    0.9343
    0.2816    0.2815    0.5657
    0.7815    0.2184    0.5657
    0.3207    0.1843    0.0657
    0.5631    0.0682    0.0657
    0.6161    0.2475    0.0657
    0.7184    0.2815    0.5657
    0.2185    0.2184    0.5657
    0.3839    0.2475    0.0657
    0.6793    0.1843    0.0657
    0.5000    0.1666    0.0657
    0.4368    0.0682    0.0657
    0.4368    0.0000    0.9343
    0.5632    0.0000    0.9343
    0.4368    0.0000    0.5657
    0.5632    0.0000    0.5657
    0.1023    0.6326    0.9343
    0.8977    0.6326    0.9343
    0.0000    0.7348    0.9343
    0.8977    0.6326    0.5657
    0.1023    0.6326    0.5657
    0.0000    0.7348    0.5657
    0.0000    0.6666    0.4343
    0.8207    0.6843    0.4343
    0.0631    0.5682    0.4343
    0.1161    0.7475    0.4343
    0.7816    0.7815    0.9343
    0.2815    0.7184    0.9343
    0.1793    0.6843    0.4343
    0.8839    0.7475    0.4343
    0.9368    0.5682    0.4343
    0.7185    0.7184    0.9343
    0.2184    0.7815    0.9343
    0.7816    0.7815    0.5657
    0.2815    0.7184    0.5657
    0.8207    0.6843    0.0657
    0.0631    0.5682    0.0657
    0.1161    0.7475    0.0657
    0.2184    0.7815    0.5657
    0.7185    0.7184    0.5657
    0.8839    0.7475    0.0657
    0.1793    0.6843    0.0657
    0.0000    0.6666    0.0657
    0.9368    0.5682    0.0657
    0.9368    0.5000    0.9343
    0.0632    0.5000    0.9343
    0.9368    0.5000    0.5657
    0.0632    0.5000    0.5657
    0.6023    0.8674    0.9343
    0.3977    0.8674    0.9343
    0.5000    0.7652    0.9343
    0.3977    0.8674    0.5657
    0.6023    0.8674    0.5657
    0.5000    0.7652    0.5657
    0.5000    0.8334    0.4343
    0.3207    0.8157    0.4343
    0.5631    0.9318    0.4343
    0.6161    0.7525    0.4343
    0.6793    0.8157    0.4343
    0.3839    0.7525    0.4343
    0.4368    0.9318    0.4343
    0.3207    0.8157    0.0657
    0.5631    0.9318    0.0657
    0.6161    0.7525    0.0657
    0.3839    0.7525    0.0657
    0.6793    0.8157    0.0657
    0.5000    0.8334    0.0657
    0.4368    0.9318    0.0657
    0.1023    0.3674    0.9343
    0.8977    0.3674    0.9343
    0.0000    0.2652    0.9343
    0.8977    0.3674    0.5657
    0.1023    0.3674    0.5657
    0.0000    0.2652    0.5657
    0.0000    0.3334    0.4343
    0.8207    0.3157    0.4343
    0.0631    0.4318    0.4343
    0.1161    0.2525    0.4343
    0.1793    0.3157    0.4343
    0.8839    0.2525    0.4343
    0.9368    0.4318    0.4343
    0.8207    0.3157    0.0657
    0.0631    0.4318    0.0657
    0.1161    0.2525    0.0657
    0.8839    0.2525    0.0657
    0.1793    0.3157    0.0657
    0.0000    0.3334    0.0657
    0.9368    0.4318    0.0657
"""
