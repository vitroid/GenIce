from dataclasses import dataclass
from typing import Optional, List, Dict, Set
from collections import defaultdict
from logging import getLogger

import numpy as np

from genice2.cell import Cell
from genice2.decorators import banner, timeit


def replicate_positions(positions1, replica_vectors, grand_cellmat):
    # レプリカ単位胞の数だけ、水分子位置を複製する。
    inv = np.linalg.inv(grand_cellmat)
    reppositions = []
    for replica_vector in replica_vectors:
        # レプリカ単位胞内での分子の位置
        replica = positions1 + replica_vector
        # 大セル内の位置に換算する
        replica = replica @ inv
        replica -= np.floor(replica)
        # 束ねて
        reppositions.append(replica)
    # 一つの配列にまとめ
    reppositions = np.vstack(reppositions)
    return reppositions


@dataclass
class Stage1Output:
    """Stage1の出力データ"""

    reppositions: np.ndarray  # 複製された分子位置
    repcell: Cell  # 複製されたセル
    repcagetype: Optional[List[str]]  # 複製されたケージタイプ
    repcagepos: Optional[np.ndarray]  # 複製されたケージ位置
    cagetypes: Optional[Dict[str, Set[int]]]  # ケージタイプの集合


class Stage1:
    """単位胞の複製を行うステージ"""

    def __init__(
        self,
        waters1: np.ndarray,
        replica_vectors: np.ndarray,
        reshape_matrix: np.ndarray,
        cell1: Cell,
        cagepos1: Optional[np.ndarray] = None,
        cagetype1: Optional[List[str]] = None,
    ):
        self.waters1 = waters1
        self.replica_vectors = replica_vectors
        self.reshape_matrix = reshape_matrix
        self.cell1 = cell1
        self.cagepos1 = cagepos1
        self.cagetype1 = cagetype1

    @timeit
    @banner
    def execute(self, noise: float = 0.0, assess_cages: bool = False) -> Stage1Output:
        """Replicates the unit cell."""
        logger = getLogger()
        reppositions = replicate_positions(
            self.waters1, self.replica_vectors, self.reshape_matrix
        )
        repcell = Cell(self.reshape_matrix @ self.cell1.mat)
        if noise > 0.0:
            logger.info(f"  Add noise: {noise}.")
            perturb = np.random.normal(
                loc=0.0,
                scale=noise * 0.01 * 3.0 * 0.5,  # in percent, radius of water
                size=reppositions.shape,
            )
            reppositions += repcell.abs2rel(perturb)

        if assess_cages:
            logger.info("  Assessing the cages...")
            self.cagepos1, self.cagetype1 = cage.assess_cages(self.graph1, self.waters1)
            logger.info("  Done assessment.")

        if self.cagepos1 is not None:
            repcagepos = replicate_positions(
                self.cagepos1, self.replica_vectors, self.reshape_matrix
            )
            repcagetype = [
                self.cagetype1[i % len(self.cagetype1)] for i in range(len(repcagepos))
            ]
            cagetypes = defaultdict(set)
            for i, typ in enumerate(repcagetype):
                cagetypes[typ].add(i)
        else:
            repcagepos = None
            repcagetype = None
            cagetypes = None

        return Stage1Output(
            reppositions=reppositions,
            repcell=repcell,
            repcagetype=repcagetype,
            repcagepos=repcagepos,
            cagetypes=cagetypes,
        )
