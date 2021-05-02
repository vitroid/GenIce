"""

Command line: /Volumes/workarea/venvs/genice2/bin/genice zeolite[ACO] -f python
Reshaping the unit cell.
  i:[1 0 0]
  j:[0 1 0]
  k:[0 0 1]
"""

import genice2.lattices
desc = {
    "ref": {
        "engel03": "Engel 2018",
        "ACO": "IZA Database"
    },
    "usage": "No options available.",
    "brief": "Hypothetical zeolitic ice"
}


class Lattice(genice2.lattices.Lattice):
    def __init__(self):
        self.bondlen = 0.30360000000000015
        self.coord = 'relative'
        from genice2.cell import cellvectors
        self.cell = cellvectors(a=0.88235294, b=0.88235294, c=0.88235294)
        self.density = 0.6961850990811413
        self.waters = """
    0.3436    0.3436    0.6564
    0.8436    0.8436    0.1564
    0.3436    0.3436    0.3436
    0.8436    0.8436    0.8436
    0.6564    0.3436    0.3436
    0.1564    0.8436    0.8436
    0.6564    0.3436    0.6564
    0.1564    0.8436    0.1564
    0.3436    0.6564    0.3436
    0.8436    0.1564    0.8436
    0.3436    0.6564    0.6564
    0.8436    0.1564    0.1564
    0.6564    0.6564    0.3436
    0.1564    0.1564    0.8436
    0.6564    0.6564    0.6564
    0.1564    0.1564    0.1564
"""
