"""
GenIce class
"""

import itertools as it
from collections import defaultdict
from logging import getLogger
from typing import Type, Union, List, Dict, Set, Optional, Tuple
from dataclasses import dataclass, field
import sys

import numpy as np
import pairlist as pl
import networkx as nx

# from genice2 import cage

# from genice2 import digraph_unused as dg
from genice2.cell import Cell, rel_wrap, cellshape
from genice2.lattices import Lattice
import genice_core

# A virtual monatomic molecule
from genice2.molecules import Molecule
from genice2.valueparser import (
    parse_cages,
    parse_pairs,
    # plugin_option_parser,
    put_in_array,
)
from genice2.stage1 import Stage1
from genice2.stage2 import Stage2, grandcell_wrap
from genice2.stage4 import Stage4
from genice2.stage5 import Stage5
from genice2.stage6 import Stage6
from genice2.stage7 import Stage7
from genice2 import ConfigurationError


# def ice_rule(dg: nx.DiGraph, strict=False) -> bool:
#     if strict:
#         for node in dg:
#             if dg.in_degree(node) != 2 or dg.out_degree(node) != 2:
#                 return False
#     else:
#         for node in dg:
#             if dg.in_degree(node) > 2 or dg.out_degree(node) > 2:
#                 return False
#     return True


def shortest_distance(coord, cell, pairs=None):
    dmin = 1e99

    if pairs is None:
        iter = it.combinations(coord, 2)
    else:
        iter = [(coord[i], coord[j]) for i, j in pairs]

    for c1, c2 in iter:
        r = cell.rel2abs(rel_wrap(c1 - c2))
        rr = r @ r

        if rr < dmin:
            dmin = rr

    return dmin**0.5


@dataclass
class GenIceConfig:
    """GenIceの設定を保持するデータクラス"""

    water: Molecule
    guests: Dict[str, Dict[str, float]] = field(default_factory=dict)
    density: float = 0
    rep: Optional[Tuple[int, int, int]] = None
    reshape: np.ndarray = field(default_factory=lambda: np.eye(3, dtype=int))
    cations: Dict[int, str] = field(default_factory=dict)
    anions: Dict[int, str] = field(default_factory=dict)
    spot_guests: Dict[int, str] = field(default_factory=dict)
    spot_groups: Dict[int, str] = field(default_factory=dict)
    asis: bool = False
    shift: Tuple[float, float, float] = (0.0, 0.0, 0.0)
    depol: str = "strict"
    assess_cages: bool = False
    noise: float = 0.0


@dataclass
class GenIceState:
    """GenIceの状態を保持するデータクラス"""

    waters1: np.ndarray
    cell1: Cell
    pairs1: Optional[List[Tuple[int, int]]] = None
    bondlen: Optional[float] = None
    cagepos1: Optional[np.ndarray] = None
    cagetype1: Optional[List[str]] = None
    fixed1: List[Tuple[int, int]] = field(default_factory=list)
    dopants1: Set[int] = field(default_factory=set)
    # groups1: Dict[int, Dict[int, str]] = field(
    #     default_factory=lambda: defaultdict(dict)
    # )
    filled_cages: Set[int] = field(default_factory=set)
    graph1: nx.Graph = field(default_factory=nx.Graph)


class GenIce:
    """氷の結晶構造を生成するクラス"""

    def __init__(
        self,
        lat: Type[Lattice],
        config: Optional[GenIceConfig] = None,
    ):
        logger = getLogger()
        logger.info(sys.argv)
        self.config = config or GenIceConfig()
        self.state = self._initialize_state(lat)
        self._setup_replica_vectors()
        self._setup_waters(lat)
        self._setup_bond_length(lat)
        self._setup_pairs(lat)
        self._setup_density(lat)
        self._setup_cages(lat)
        self._setup_fixed_bonds(lat)
        # self._setup_dopants(lat)
        self._setup_graph(lat)
        self._setup_groups(lat)

    def _initialize_state(self, lat: Type[Lattice]) -> GenIceState:
        """初期状態を設定"""
        return GenIceState(
            waters1=np.array([]),
            cell1=Cell(lat.cell),
        )

    def _setup_replica_vectors(self):
        """レプリカベクトルの設定"""
        logger = getLogger()
        if self.config.rep is not None:
            logger.warning("rep for GenIce() is deprecated. Use reshape instead.")
            self.replica_vectors = np.array(
                [
                    (x, y, z)
                    for x in range(self.config.rep[0])
                    for y in range(self.config.rep[1])
                    for z in range(self.config.rep[2])
                ]
            )
            self.reshape_matrix = np.diag(self.config.rep)
        else:
            self._setup_reshape_matrix()

    def _setup_reshape_matrix(self):
        """reshape行列の設定"""
        logger = getLogger()
        logger.info("  Reshaping the unit cell.")
        self.reshape_matrix = self.config.reshape

        i, j, k = np.array(self.reshape_matrix)
        logger.info(f"    i:{i}")
        logger.info(f"    j:{j}")
        logger.info(f"    k:{k}")

        a, b, c, A, B, C = cellshape(self.reshape_matrix @ self.state.cell1.mat)
        logger.info("  Reshaped cell:")
        logger.info(f"    a,b,c = {a}, {b}, {c}")
        logger.info(f"    A,B,C = {A}, {B}, {C}")

        self.replica_vectors = self._calculate_replica_vectors(i, j, k)

        # レプリカ単位胞には0から順に番号ラベルがついている。replica_vector_labelsは単位胞の位置をラベルに変換する
        self.replica_vector_labels = {
            tuple(xyz): i for i, xyz in enumerate(self.replica_vectors)
        }

    def _calculate_replica_vectors(
        self, i: np.ndarray, j: np.ndarray, k: np.ndarray
    ) -> np.ndarray:
        """レプリカベクトルの計算"""
        logger = getLogger()
        corners = np.array(
            [a * i + b * j + c * k for a in (0, 1) for b in (0, 1) for c in (0, 1)]
        )

        mins = np.min(corners, axis=0)
        maxs = np.max(corners, axis=0)
        logger.info(f"mins: {mins}, maxs: {maxs}")

        det = abs(np.linalg.det(self.reshape_matrix))
        det = np.floor(det + 0.5).astype(int)
        invdet = np.floor(np.linalg.inv(self.reshape_matrix) * det + 0.5).astype(int)

        vecs = set()
        for a in range(mins[0], maxs[0] + 1):
            for b in range(mins[1], maxs[1] + 1):
                for c in range(mins[2], maxs[2] + 1):
                    abc = np.array([a, b, c])
                    rep = grandcell_wrap(abc, self.reshape_matrix, invdet, det).astype(
                        int
                    )
                    if tuple(rep) not in vecs:
                        vecs.add(tuple(rep))

        vecs = np.array(list(vecs))
        vol = abs(np.linalg.det(self.reshape_matrix))
        assert np.allclose(vol, len(vecs)), (vol, vecs)

        return vecs

    def _setup_waters(self, lat: Type[Lattice]):
        """水分子の設定"""
        logger = getLogger()
        self.state.waters1 = put_in_array(lat.waters)
        logger.debug(f"Waters: {len(self.state.waters1)}")
        self.state.waters1 = self.state.waters1.reshape((-1, 3))

        if lat.coord == "absolute":
            self.state.waters1 = self.state.cell1.abs2rel(self.state.waters1)

        self.state.waters1 = np.array(self.state.waters1) + np.array(self.config.shift)
        self.state.waters1 -= np.floor(self.state.waters1)

    def _setup_pairs(self, lat: Type[Lattice]):
        """分子間結合の設定"""
        logger = getLogger()

        if "pairs" in lat.__dict__ and lat.pairs is not None:
            self.state.pairs1 = parse_pairs(lat.pairs)
        else:
            logger.info("  Pairs are not given explicitly.")
            logger.info("  Estimating the bonds according to the pair distances.")

            logger.debug(f"Bondlen: {self.state.bondlen}")
            # make bonded pairs according to the pair distance.
            # make before replicating them.
            self.state.pairs1 = [
                (i, j)
                for i, j in pl.pairs_iter(
                    self.state.waters1,
                    self.state.bondlen,
                    self.state.cell1.mat,
                    distance=False,
                )
            ]

    def _setup_bond_length(self, lat: Type[Lattice]):
        """結合長の設定"""
        logger = getLogger()
        nmol = self.state.waters1.shape[0]
        volume = self.state.cell1.volume()

        try:
            self.state.bondlen = lat.bondlen
            logger.info(f"Bond length (specified): {self.state.bondlen}")
        except AttributeError:
            logger.debug("  Estimating the bond threshold length...")
            rc = (volume / nmol) ** (1 / 3) * 1.5
            p = pl.pairs_iter(
                self.state.waters1,
                maxdist=rc,
                cell=self.state.cell1.mat,
                distance=False,
            )
            self.state.bondlen = 1.1 * shortest_distance(
                self.state.waters1, self.state.cell1, pairs=p
            )
            logger.info(f"Bond length (estim.): {self.state.bondlen}")

    def _setup_density(self, lat: Type[Lattice]):
        """密度の設定"""
        logger = getLogger()
        nmol = self.state.waters1.shape[0]
        volume = self.state.cell1.volume()
        mass = 18  # water
        NB = 6.022e23
        density0 = mass * nmol / (NB * volume * 1e-21)

        if self.config.density <= 0:
            try:
                self.density = lat.density
            except AttributeError:
                logger.info(
                    "Density is not specified. Assume the density from lattice."
                )
                dmin = shortest_distance(self.state.waters1, self.state.cell1)
                logger.info(
                    f"Closest pair distance: {dmin} (should be around 0.276 nm)"
                )
                self.density = density0 / (0.276 / dmin) ** 3
        else:
            self.density = self.config.density

        logger.info(f"Target Density: {self.density}")
        logger.info(f"Original Density: {density0}")

        ratio = (density0 / self.density) ** (1.0 / 3.0)
        self.state.cell1.scale(ratio)

        if self.state.bondlen is not None:
            self.state.bondlen *= ratio
        logger.info(f"Bond length (scaled, nm): {self.state.bondlen}")

    def _setup_cages(self, lat: Type[Lattice]):
        """ケージの設定"""
        logger = getLogger()
        if "cages" in lat.__dict__:
            self.state.cagepos1, self.state.cagetype1 = parse_cages(lat.cages)
            logger.warn("Use of 'cages' in a lattice-plugin is deprecated.")
        elif "cagepos" in lat.__dict__:
            self.state.cagepos1, self.state.cagetype1 = (
                np.array(lat.cagepos),
                lat.cagetype,
            )

        if self.state.cagepos1 is not None:
            self.state.cagepos1 = np.array(self.state.cagepos1) + np.array(
                self.config.shift
            )
            self.state.cagepos1 -= np.floor(self.state.cagepos1)

    def _setup_fixed_bonds(self, lat: Type[Lattice]):
        """固定結合の設定"""
        logger = getLogger()
        if "fixed" in lat.__dict__ and lat.fixed is not None:
            self.state.fixed1 = parse_pairs(lat.fixed)
            logger.info("Orientations of some edges are fixed.")

        if self.config.asis and len(self.state.fixed1) == 0:
            self.state.fixed1 = self.state.pairs1

    # def __setup_dopants(self, lat: Type[Lattice]):
    #     """Hook関数の設定"""
    #     if "dopeIonsToUnitCell" in lat.__dict__:
    #         self.doping_hook_function = lat.dopeIonsToUnitCell
    #     else:
    #         self.doping_hook_function = None

    def _setup_graph(self, lat: Type[Lattice]):
        """グラフの設定"""
        logger = getLogger()
        logger.info(f"Generating the graph...")

        self.state.graph1 = nx.Graph(self.state.pairs1)
        logger.info(f"graph: {self.state.graph1}")

    def _setup_groups(self, lat: Type[Lattice]):
        """グループの設定"""
        logger = getLogger()
        logger.info("Generating the groups...")
        self.state.groups1 = defaultdict(dict)

    def _requires(self, stage: int):
        logger = getLogger()
        if f"stage{stage}_output" in self.state.__dict__:
            logger.debug(f"Stage {stage} is already done.")
            return
        if stage > 1:
            self._requires(stage - 1)
        if stage == 1:
            self.state.stage1_output = Stage1(
                self.state.waters1,
                self.replica_vectors,
                self.reshape_matrix,
                self.state.cell1,
                self.state.cagepos1,
                self.state.cagetype1,
            ).execute(self.config.noise, self.config.assess_cages)
        elif stage == 2:
            self.state.stage2_output = Stage2(
                self.state.graph1,
                self.state.waters1,
                self.replica_vectors,
                self.state.fixed1,
                self.replica_vector_labels,
                self.reshape_matrix,
                self.config.anions,
                self.config.cations,
                self.state.groups1,
                self.state.cagepos1,
            ).execute()
        elif stage == 3:
            pass
        elif stage == 4:
            self.state.stage4_output = Stage4(
                self.state.stage2_output.graph,
                self.state.stage1_output.reppositions,
                self.state.stage2_output.fixed_edges,
            ).execute(self.config.depol)
        elif stage == 5:
            self.state.stage5_output = Stage5(
                self.state.stage1_output.reppositions,
                self.state.stage4_output.digraph,
                self.state.stage1_output.repcell,
                self.state.stage2_output.dopants,
            ).execute()
        elif stage == 6:
            self.state.stage6_output = Stage6(
                self.state.stage1_output.reppositions,
                self.state.stage1_output.repcell,
                self.state.stage5_output.rotmatrices,
                self.config.water,
                self.state.stage2_output.dopants,
            ).execute()
        elif stage == 7:
            self.state.stage7_output = Stage7(
                self.state.stage6_output.universe,
                self.state.stage1_output.repcagepos,
                self.state.stage1_output.repcagetype,
                self.state.stage1_output.cagetypes,
                self.state.stage2_output.filled_cages,
                self.config.guests,
                self.config.spot_guests,
                self.config.spot_groups,
                self.state.stage2_output.groups,
                self.state.stage2_output.dopants,
                self.state.stage1_output.reppositions,
                self.state.stage5_output.rotmatrices,
                self.state.stage1_output.repcell,
            ).execute()
        else:
            raise ValueError(f"Invalid stage: {stage}")

    def full_atomic_positions(self):
        """Atomic positions of all atoms."""
        self._requires(7)
        return self.state.stage7_output.universe

    def water_atomic_positions(self):
        """Atomic positions of water molecules."""
        self._requires(6)
        return self.state.stage6_output.universe

    def rotation_matrices(self):
        """Rotation matrices."""
        self._requires(5)
        return self.state.stage5_output.rotmatrices

    def cell_matrix(self):
        """Cell matrix."""
        self._requires(1)
        return self.state.stage1_output.repcell.mat

    def water_positions(self):
        """Molecular positions of water molecules in the fractional coordinates."""
        self._requires(1)
        return self.state.stage1_output.reppositions

    def hydrogen_bond_digraph(self):
        """Hydrogen bond graph."""
        self._requires(4)
        return self.state.stage4_output.digraph

    def hydrogen_bond_graph(self):
        """Hydrogen bond graph."""
        self._requires(2)
        return self.state.stage2_output.graph


def groups_info(groups):
    logger = getLogger()
    for root, cages in groups.items():
        for cage, group in cages.items():
            logger.info(f"    Group {group} of dopant {root} in cage {cage}")


def guests_info(cagetypes, molecules):
    logger = getLogger()
    for cagetype, cageid in cagetypes.items():
        logger.info(f"    Guests in cage type {cagetype}:")

        for molec, cages in molecules.items():
            cages = set(cages)
            cages &= cageid

            if len(cages):
                logger.info(f"      {molec} * {len(cages)} @ {cages}")
