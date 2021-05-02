"""

Command line: /Volumes/workarea/venvs/genice2/bin/genice zeolite[SGT] -f python
Reshaping the unit cell.
  i:[1 0 0]
  j:[0 1 0]
  k:[0 0 1]
"""
import genice2.lattices
desc = {
    "ref": {
        "engel31": "Engel 2018",
        "SGT": "IZA Database"
    },
    "usage": "No options available.",
    "brief": "Hypothetical zeolitic ice"
}


class Lattice(genice2.lattices.Lattice):
    def __init__(self):
        self.bondlen = 0.30360000000000026
        self.coord = 'relative'
        from genice2.cell import cellvectors
        self.cell = cellvectors(a=0.95247725, b=0.95247725, c=3.16883803)
        self.density = 0.6654301389641005
        self.waters = """
    0.5000    0.9031    0.2162
    0.0000    0.4031    0.7162
    0.1531    0.7500    0.9662
    0.6531    0.2500    0.4662
    0.0000    0.0969    0.7162
    0.5000    0.5969    0.2162
    0.3469    0.2500    0.4662
    0.8469    0.7500    0.9662
    0.5000    0.0969    0.7838
    0.0000    0.5969    0.2838
    0.8469    0.2500    0.0338
    0.3469    0.7500    0.5338
    0.0000    0.9031    0.2838
    0.5000    0.4031    0.7838
    0.6531    0.7500    0.5338
    0.1531    0.2500    0.0338
    0.7500    0.9676    0.7500
    0.2500    0.4676    0.2500
    0.2824    0.5000    0.0000
    0.7824    0.0000    0.5000
    0.7500    0.0324    0.2500
    0.2500    0.5324    0.7500
    0.2176    0.5000    0.5000
    0.7176    0.0000    0.0000
    0.2500    0.9676    0.7500
    0.7500    0.4676    0.2500
    0.7176    0.5000    0.0000
    0.2176    0.0000    0.5000
    0.2500    0.0324    0.2500
    0.7500    0.5324    0.7500
    0.7824    0.5000    0.5000
    0.2824    0.0000    0.0000
    0.5000    0.9012    0.0547
    0.0000    0.4012    0.5547
    0.1512    0.7500    0.8047
    0.6512    0.2500    0.3047
    0.0000    0.0988    0.5547
    0.5000    0.5988    0.0547
    0.3488    0.2500    0.3047
    0.8488    0.7500    0.8047
    0.5000    0.0988    0.9453
    0.0000    0.5988    0.4453
    0.8488    0.2500    0.1953
    0.3488    0.7500    0.6953
    0.0000    0.9012    0.4453
    0.5000    0.4012    0.9453
    0.6512    0.7500    0.6953
    0.1512    0.2500    0.1953
    0.5000    0.0360    0.1318
    0.0000    0.5360    0.6318
    0.2860    0.7500    0.8818
    0.7860    0.2500    0.3818
    0.0000    0.9640    0.6318
    0.5000    0.4640    0.1318
    0.2140    0.2500    0.3818
    0.7140    0.7500    0.8818
    0.5000    0.9640    0.8682
    0.0000    0.4640    0.3682
    0.7140    0.2500    0.1182
    0.2140    0.7500    0.6182
    0.0000    0.0360    0.3682
    0.5000    0.5360    0.8682
    0.7860    0.7500    0.6182
    0.2860    0.2500    0.1182
"""
