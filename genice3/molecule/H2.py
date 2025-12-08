import genice3.molecule
import numpy as np

desc = {
    "ref": {"H2": "https://www.britannica.com/science/hydrogen"},
    "brief": "Hydrogen molecule.",
    "usage": "No options available.",
}


class Molecule(genice3.molecule.Molecule):
    """
    水素分子を定義するクラス。
    """

    def __init__(self):
        sites = np.array([[0, 0, -0.037], [0, 0, +0.037]])  # nm, HH
        labels = ["H", "H"]
        name = "H2"
        is_water = False
        super().__init__(sites=sites, labels=labels, name=name, is_water=is_water)

