import genice3.molecule
import numpy as np
import math

desc = {
    "ref": {
        "TIP3P": "Jorgensen, Chandrasekhar, Madura, Impey, Klein, J Chem Phys, 79, 926 (1983)."
    },
    "brief": "Physical model of water; Oxygen atom is on the lattice point.",
    "usage": "No options available.",
}


class Molecule(genice3.molecule.Molecule):
    """
    物理的な水分子モデルを定義するクラス。酸素原子が格子点上にある。
    """

    def __init__(self):
        oh = 0.09572
        hangle = 104.52 * math.pi / 180 / 2
        mass = 18
        ohz = oh * math.cos(hangle)
        ohy = oh * math.sin(hangle)
        oz = 0  # -ohz*2/mass
        sites = np.array(
            [[0, 0, oz], [0, ohy, ohz + oz], [0, -ohy, ohz + oz]]
        )  # nm, OHHM

        labels = ["O", "H", "H"]
        name = "SOL"
        is_water = True
        super().__init__(sites=sites, labels=labels, name=name, is_water=is_water)
