"""

Command line: /Volumes/workarea/venvs/genice2/bin/genice zeolite[MAR] -f python
Reshaping the unit cell.
  i:[1 0 0]
  j:[0 1 0]
  k:[0 0 1]
"""
import genice2.lattices
desc = {
    "ref": {
        "engel29": "Engel 2018",
        "MAR": "IZA Database"
    },
    "usage": "No options available.",
    "brief": "Hypothetical zeolitic ice"
}


class Lattice(genice2.lattices.Lattice):
    def __init__(self):
        self.bondlen = 0.30360000000000026
        self.coord = 'relative'

        import numpy as np
        self.cell = np.array([[1.16251300,
                               0.00000000,
                               0.00000000],
                              [-0.58125650,
                               1.00676579,
                               0.00000000],
                              [0.00000000,
                               0.00000000,
                               2.85834478],
                              ])
        self.density = 0.6433145477702799
        self.waters = """
    0.7506    0.7506    0.5000
    0.2494    0.0000    0.5000
    1.0000    0.2494    0.5000
    0.2494    0.2494    0.0000
    0.7506    0.0000    0.0000
    0.0000    0.7506    0.0000
    0.2494    0.2494    0.5000
    0.0000    0.7506    0.5000
    0.7506    0.0000    0.5000
    0.7506    0.7506    0.0000
    1.0000    0.2494    0.0000
    0.2494    0.0000    0.0000
    0.0837    0.6701    0.7500
    0.3299    0.4136    0.7500
    0.5864    0.9163    0.7500
    0.9163    0.3299    0.2500
    0.6701    0.5864    0.2500
    0.4136    0.0837    0.2500
    0.3299    0.9163    0.7500
    0.5864    0.6701    0.7500
    0.0837    0.4136    0.7500
    0.6701    0.0837    0.2500
    0.4136    0.3299    0.2500
    0.9163    0.5864    0.2500
    0.5876    0.9229    0.3282
    0.0771    0.6647    0.3282
    0.3353    0.4124    0.3282
    0.4124    0.0771    0.8282
    0.9229    0.3353    0.8282
    0.6647    0.5876    0.8282
    0.0771    0.4124    0.3282
    0.3353    0.9229    0.3282
    0.5876    0.6647    0.3282
    0.9229    0.5876    0.8282
    0.6647    0.0771    0.8282
    0.4124    0.3353    0.8282
    0.4124    0.0771    0.6718
    0.9229    0.3353    0.6718
    0.6647    0.5876    0.6718
    0.5876    0.9229    0.1718
    0.0771    0.6647    0.1718
    0.3353    0.4124    0.1718
    0.9229    0.5876    0.6718
    0.6647    0.0771    0.6718
    0.4124    0.3353    0.6718
    0.0771    0.4124    0.1718
    0.3353    0.9229    0.1718
    0.5876    0.6647    0.1718
    0.9156    0.3319    0.4131
    0.6681    0.5837    0.4131
    0.4163    0.0844    0.4131
    0.0844    0.6681    0.9131
    0.3319    0.4163    0.9131
    0.5837    0.9156    0.9131
    0.6681    0.0844    0.4131
    0.4163    0.3319    0.4131
    0.9156    0.5837    0.4131
    0.3319    0.9156    0.9131
    0.5837    0.6681    0.9131
    0.0844    0.4163    0.9131
    0.0844    0.6681    0.5869
    0.3319    0.4163    0.5869
    0.5837    0.9156    0.5869
    0.9156    0.3319    0.0869
    0.6681    0.5837    0.0869
    0.4163    0.0844    0.0869
    0.3319    0.9156    0.5869
    0.5837    0.6681    0.5869
    0.0844    0.4163    0.5869
    0.6681    0.0844    0.0869
    0.4163    0.3319    0.0869
    0.9156    0.5837    0.0869
"""