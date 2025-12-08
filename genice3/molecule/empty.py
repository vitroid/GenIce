import numpy as np
import genice3.molecule


class Molecule(genice3.molecule.Molecule):
    """
    空のケージ用のプレースホルダー。
    """

    def __init__(self):
        # placeholder for empty cage

        sites = []
        labels = []
        name = "empty"
        is_water = False
        super().__init__(sites=sites, labels=labels, name=name, is_water=is_water)
