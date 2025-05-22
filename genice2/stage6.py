from dataclasses import dataclass
from typing import List, Dict

import numpy as np

from genice2.cell import Cell
from genice2.decorators import banner, timeit
from genice2.molecules import Molecule, AtomicStructure


@dataclass
class Stage6Output:
    """Stage6の出力データ"""

    universe: List[np.ndarray]  # 原子位置


class Stage6:
    """水分子と置換分子の原子配置を行うステージ"""

    def __init__(
        self,
        reppositions: np.ndarray,
        repcell: Cell,
        rotmatrices: np.ndarray,
        water: Molecule,
        dopants: Dict[int, str],
    ):
        self.reppositions = reppositions
        self.repcell = repcell
        self.rotmatrices = rotmatrices
        self.water = water
        self.dopants = dopants

    @timeit
    @banner
    def execute(self) -> Stage6Output:
        """Arrange atoms of water and replacements"""
        universe = []
        universe.append(
            AtomicStructure(
                self.reppositions,
                self.repcell,
                self.rotmatrices,
                self.water,
                immutables=set(self.dopants),
            )
        )
        return Stage6Output(universe=universe)
