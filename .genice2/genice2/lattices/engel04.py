"""

Command line: /Volumes/workarea/venvs/genice2/bin/genice zeolite[LTA] -f python
Reshaping the unit cell.
  i:[1 0 0]
  j:[0 1 0]
  k:[0 0 1]
"""
import genice2.lattices
desc = {
    "ref": {
        "engel04": "Engel 2018",
        "LTA": "IZA Database"
    },
    "usage": "No options available.",
    "brief": "Hypothetical zeolitic ice"
}


class Lattice(genice2.lattices.Lattice):
    def __init__(self):
        self.bondlen = 0.30360000000000004
        self.coord = 'relative'
        from genice2.cell import cellvectors
        self.cell = cellvectors(a=1.07055113, b=1.07055113, c=1.07055113)
        self.density = 0.5846833799862865
        self.waters = """
    0.6316    0.1823    0.0000
    0.8177    0.0000    0.3684
    0.6316    0.8177    0.0000
    0.8177    0.0000    0.6316
    0.1823    0.6316    0.0000
    0.0000    0.8177    0.3684
    0.8177    0.6316    0.0000
    0.0000    0.8177    0.6316
    0.0000    0.6316    0.8177
    0.3684    0.8177    0.0000
    0.0000    0.6316    0.1823
    0.0000    0.3684    0.8177
    0.1823    0.0000    0.6316
    0.3684    0.1823    0.0000
    0.0000    0.3684    0.1823
    0.1823    0.0000    0.3684
    0.6316    0.0000    0.8177
    0.8177    0.3684    0.0000
    0.6316    0.0000    0.1823
    0.3684    0.0000    0.8177
    0.0000    0.1823    0.6316
    0.1823    0.3684    0.0000
    0.3684    0.0000    0.1823
    0.0000    0.1823    0.3684
"""
