import genice3.molecule
import numpy as np
import math

desc = {
    "ref": {"TIP7P": "Zhao 2019"},
    "brief": "A seven-site water model.",
    "usage": "No options available.",
}


class Molecule(genice3.molecule.Molecule):
    """
    7サイトモデルの水分子を定義するクラス。
    """

    def __init__(self):
        oh = 0.04786 * 2  # nm
        ol = 0.041  # nm
        hangle = math.radians(104.52) / 2
        mangle = math.radians(109.47) / 2
        mass = 18
        ohz = oh * math.cos(hangle)
        ohy = oh * math.sin(hangle)
        olz = -ol * math.cos(mangle)
        olx = ol * math.sin(mangle)
        oz = -ohz * 2 / mass

        sites = np.array(
            [
                [0, 0, oz],
                [0, ohy, ohz + oz],
                [0, -ohy, ohz + oz],
                [0, ohy / 2, ohz / 2 + oz],
                [0, -ohy / 2, ohz / 2 + oz],
                [olx, 0, olz + oz],
                [-olx, 0, olz + oz],
            ]
        )  # nm, OHHMMLL

        labels = ["OW", "HW1", "HW2", "MW1", "MW2", "LW1", "LW2"]
        name = "SOL"
        is_water = True
        super().__init__(sites=sites, labels=labels, name=name, is_water=is_water)

