import math
import numpy as np
import genice3.molecule

desc = {
    "usage": "No options available.",
    "brief": "A typical 3-site model.",
}


class Molecule(genice3.molecule.Molecule):
    """
    SPC/Eモデルの水分子を定義するクラス。
    """

    def __init__(self):
        oh = 0.10000
        hangle = 109.47 * math.pi / 180 / 2
        mass = 18
        ohz = oh * math.cos(hangle)
        ohy = oh * math.sin(hangle)
        oz = -ohz * 2 / mass

        sites = np.array([[0, 0, oz], [0, ohy, ohz + oz], [0, -ohy, ohz + oz]])  # nm

        labels = ["Ow", "Hw", "Hw"]
        name = "SOL"
        is_water = True
        super().__init__(sites=sites, labels=labels, name=name, is_water=is_water)
