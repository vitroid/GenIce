"""

Command line: /Volumes/workarea/venvs/genice2/bin/genice zeolite[LTA] -f python
Reshaping the unit cell.
  i:[1 0 0]
  j:[0 1 0]
  k:[0 0 1]
"""
desc={
    "ref": {
        "engel04": "Engel, E.A., Anelli, A., Ceriotti, M. et al. Mapping uncharted territory in ice from zeolite networks to ice structures. Nat Commun 9, 2173 (2018). https://doi.org/10.1038/s41467-018-04618-6",
        "LTA": "Database of Zeolite Structures, https://asia.iza-structure.org/IZA-SC/framework.php?STC=LTA"
    },
    "usage": "No options available.",
    "brief": "Hypothetical zeolitic ice"
}

import genice2.lattices
class Lattice(genice2.lattices.Lattice):
    def __init__(self):
        self.bondlen=0.30360000000000004
        self.coord='relative'
        from genice2.cell import cellvectors
        self.cell = cellvectors(a=1.07055113, b=1.07055113, c=1.07055113)
        self.density=0.5846833799862865
        self.waters="""
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
