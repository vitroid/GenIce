
desc={"ref": {},
      "usage": "No options available.",
      "brief": "Half lattice of ice VI."
      }


bondlen=2.3681227356441177
coord='relative'
density=1.373/2
waters="""
    0.2200    0.5000    0.3800
    0.7800    0.5000    0.3800
    0.5000    0.2200    0.6200
    0.5000    0.5000    0.0000
    0.5000    0.7800    0.6200
"""

cages="""
Oc    0.0000    0.0000    0.0000
"""


from genice.cell import cellvectors
cell = cellvectors(a=4.87672629,
                   b=4.87385128,
                   c=4.49131038)
