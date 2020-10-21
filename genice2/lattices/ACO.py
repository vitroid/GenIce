"""

Command line: /Volumes/workarea/venvs/genice2/bin/genice zeolite[ACO] -f python
Reshaping the unit cell.
  i:[1 0 0]
  j:[0 1 0]
  k:[0 0 1]
"""

desc={
    "ref": {
        "engel03": "Engel, E.A., Anelli, A., Ceriotti, M. et al. Mapping uncharted territory in ice from zeolite networks to ice structures. Nat Commun 9, 2173 (2018). https://doi.org/10.1038/s41467-018-04618-6",
        "ACO": "Database of Zeolite Structures, https://asia.iza-structure.org/IZA-SC/framework.php?STC=ACO"
    },
    "usage": "No options available.",
    "brief": "Hypothetical zeolitic ice"
}
import genice2.lattices
class Lattice(genice2.lattices.Lattice):
    def __init__(self):
        self.bondlen=0.30360000000000015
        self.coord='relative'
        from genice2.cell import cellvectors
        self.cell = cellvectors(a=0.88235294, b=0.88235294, c=0.88235294)
        self.density=0.6961850990811413
        self.waters="""
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
