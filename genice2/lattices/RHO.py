# coding: utf-8
"""
Data source:
Huang, Y et al. “A New Phase Diagram of Water Under Negative Pressure: the Rise of the Lowest-Density Clathrate S-III.” Science Advances 2.2 (2016): e1501010–e1501010.
"""
from genice2.cell import cellvectors
import genice2.lattices
desc = {"ref": {"sIII": 'Huang 2016',
                "RHO": 'IZA Database'},
        "usage": "No options available.",
        "brief": "Hypothetical ice at negative pressure ice 'sIII'."
        }


class Lattice(genice2.lattices.Lattice):
    def __init__(self):
        self.cell = """
        13.339813507 13.339813507 13.339813507
        """
        self.waters = """
        8.05324541418 1.38333866068 10.0048601303
        1.38333866068 8.05324541418 3.33495337675
        11.9564748463 3.33495337675 8.05324541418
        5.28656809283 10.0048601303 1.38333866068
        10.0048601303 5.28656809283 11.9564748463
        3.33495337675 11.9564748463 5.28656809283
        8.05324541418 1.38333866068 3.33495337675
        1.38333866068 8.05324541418 10.0048601303
        11.9564748463 3.33495337675 5.28656809283
        5.28656809283 10.0048601303 11.9564748463
        10.0048601303 5.28656809283 1.38333866068
        3.33495337675 11.9564748463 8.05324541418
        11.9564748463 8.05324541418 3.33495337675
        5.28656809283 1.38333866068 10.0048601303
        10.0048601303 11.9564748463 5.28656809283
        3.33495337675 5.28656809283 11.9564748463
        8.05324541418 10.0048601303 1.38333866068
        1.38333866068 3.33495337675 8.05324541418
        11.9564748463 8.05324541418 10.0048601303
        5.28656809283 1.38333866068 3.33495337675
        10.0048601303 11.9564748463 8.05324541418
        3.33495337675 5.28656809283 1.38333866068
        8.05324541418 10.0048601303 11.9564748463
        1.38333866068 3.33495337675 5.28656809283
        5.28656809283 11.9564748463 3.33495337675
        11.9564748463 5.28656809283 10.0048601303
        1.38333866068 10.0048601303 5.28656809283
        8.05324541418 3.33495337675 11.9564748463
        3.33495337675 8.05324541418 1.38333866068
        10.0048601303 1.38333866068 8.05324541418
        5.28656809283 11.9564748463 10.0048601303
        11.9564748463 5.28656809283 3.33495337675
        1.38333866068 10.0048601303 8.05324541418
        8.05324541418 3.33495337675 1.38333866068
        3.33495337675 8.05324541418 11.9564748463
        10.0048601303 1.38333866068 5.28656809283
        8.05324541418 11.9564748463 3.33495337675
        1.38333866068 5.28656809283 10.0048601303
        11.9564748463 10.0048601303 5.28656809283
        5.28656809283 3.33495337675 11.9564748463
        10.0048601303 8.05324541418 1.38333866068
        3.33495337675 1.38333866068 8.05324541418
        8.05324541418 11.9564748463 10.0048601303
        1.38333866068 5.28656809283 3.33495337675
        11.9564748463 10.0048601303 8.05324541418
        5.28656809283 3.33495337675 1.38333866068
        10.0048601303 8.05324541418 11.9564748463
        3.33495337675 1.38333866068 5.28656809283
        """
        self.coord = "absolute"
        self.bondlen = 3
        self.density = 0.604398971981

        self.cell = cellvectors(a=13.339813507,
                                b=13.339813507,
                                c=13.339813507)
