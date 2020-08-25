# coding: utf-8
"""
Data sources

"""

desc={"ref": {"Ic": 'W. L. Vos, L. W. Finger, R. J. Hemley, H.-K. Mao, Phys. Rev. Lett., 1993, DOI:10.1103/PhysRevLett.71.3150.'},
      "usage": "No options available.",
      "brief": "Ice Ic."
      }

from genice.cell import cellvectors
import genice.lattices

class Lattice(genice.lattices.Lattice):
    def __init__(self):
        self.cell = cellvectors(a=4.0,
                           b=4.0,
                           c=4.0)

        self.density = 0.92     #default density

        self.bondlen = 1.9      #bond threshold

        self.waters = """
        0 0 0
        0.5 0.5 0
        0.5 0 0.5
        0 0.5 0.5
        0.25 0.25 0.25
        0.75 0.75 0.25
        0.75 0.25 0.75
        0.25 0.75 0.75
        """
        self.coord = "relative"

        self.pairs="""
        0 4
        0 5
        0 6
        0 7
        1 4
        1 5
        1 6
        1 7
        2 4
        2 5
        2 6
        2 7
        3 4
        3 5
        3 6
        3 7
        """
