from genice2.cell import cellvectors
import genice2.lattices

desc = {
    "ref": {},
    "usage": "No options available.",
    "brief": "Half lattice of ice VI.",
    "test": ({"options": "-r 2 2 2"},),
}


class Lattice(genice2.lattices.Lattice):
    def __init__(self):
        self.bondlen = 2.3681227356441177
        self.coord = "relative"
        self.density = 1.373 / 2
        self.waters = """
            0.2200    0.5000    0.3800
            0.7800    0.5000    0.3800
            0.5000    0.2200    0.6200
            0.5000    0.5000    0.0000
            0.5000    0.7800    0.6200
        """

        self.cages = """
        Oc    0.0000    0.0000    0.0000
        """

        self.cell = cellvectors(a=4.87672629, b=4.87385128, c=4.49131038)
