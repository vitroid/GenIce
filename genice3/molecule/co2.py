import numpy as np
import genice3.molecule

desc = {"usage": "No options available.", "brief": "An all-atom methane model."}


class Molecule(genice3.molecule.Molecule):
    """
    二酸化炭素分子を定義するクラス。
    """

    def __init__(self):
        sites = np.array([[0, 0, 0], [0, 0, -0.1149], [0, 0, +0.1149]])  # nm, OHHM

        labels = ["C", "O", "O"]
        name = "CO2"
        is_water = False
        super().__init__(sites=sites, labels=labels, name=name, is_water=is_water)
