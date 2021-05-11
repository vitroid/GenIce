# coding: utf-8

from genice2.cell import cellvectors
import genice2.lattices
desc = {"ref": {"0": "Russo 2014",
                },
        "usage": "No options available.",
        "brief": 'Metastable ice "0".'
        }

# Keep this file as simple as possible for compatibility


class Lattice(genice2.lattices.Lattice):
    def __init__(self):
        self.density = 0.92  # default self.density

        self.bondlen = 3  # bond threshold
        self.cell = """
        5.85695813012517 5.85695813012517 10.4558393727084
        """

        self.waters = """
        4.89641699678464 3.88902019840311 8.96065434241106
        3.88902019840311 4.89641699678464 6.72310471665148
        0.960541133340528 3.88902019840311 3.73273465605689
        0 0 5.22791968635418
        1.96793793172206 0.960541133340528 6.72310471665148
        2.92847906506258 2.92847906506258 0
        4.89641699678464 1.96793793172206 3.73273465605689
        1.96793793172206 4.89641699678464 1.4951850302973
        0 0 0
        3.88902019840311 0.960541133340528 1.4951850302973
        0.960541133340528 1.96793793172206 8.96065434241106
        2.92847906506258 2.92847906506258 5.22791968635418
        """

        self.coord = "absolute"

        self.cell = cellvectors(a=5.85695813012517,
                                b=5.85695813012517,
                                c=10.4558393727084)
