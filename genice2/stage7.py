from dataclasses import dataclass
from typing import List, Dict, Set
from logging import getLogger
from collections import defaultdict

import numpy as np

from genice2.decorators import banner, timeit
from genice2.molecules import AtomicStructure
from genice2.plugin import Group, safe_import
from genice2.valueparser import plugin_option_parser
from genice2.cell import Cell, rel_wrap


def guests_info(cagetypes, molecules):
    logger = getLogger()
    for cagetype, cageid in cagetypes.items():
        logger.info(f"    Guests in cage type {cagetype}:")

        for molec, cages in molecules.items():
            cages = set(cages)
            cages &= cageid

            if len(cages):
                logger.info(f"      {molec} * {len(cages)} @ {cages}")


def groups_info(groups):
    logger = getLogger()
    for root, cages in groups.items():
        for cage, group in cages.items():
            logger.info(f"    Group {group} of dopant {root} in cage {cage}")


def neighbor_cages_of_dopants(dopants, waters, cagepos, cell):
    """
    Just shows the environments of the dopants
    """
    # logger = getLogger()
    dnei = defaultdict(set)

    for site, name in dopants.items():
        org = waters[site]

        for i, pos in enumerate(cagepos):
            # Displacement (relative)
            a = cell.rel2abs(rel_wrap(pos - org))
            sqdistance = a @ a

            if sqdistance < 0.57**2:
                dnei[site].add(i)

    return dnei


def dopants_info(dopants=None, waters=None, cagepos=None, cell=None):
    logger = getLogger()
    if dopants is None:
        dopants = self.dopants

    if waters is None:
        waters = self.waters

    if cagepos is None:
        cagepos = self.cagepos

    if cell is None:
        cell = self.cell

    dopants_neighbors = neighbor_cages_of_dopants(dopants, waters, cagepos, cell)

    for dopant, cages in dopants_neighbors.items():
        logger.info(f"    Cages adjacent to dopant {dopant}: {cages}")

    return dopants_neighbors


@dataclass
class Stage7Output:
    """Stage7の出力データ"""

    universe: List[np.ndarray]  # 原子位置（ゲスト分子を含む）


class Stage7:
    """ゲスト分子の配置を行うステージ"""

    def __init__(
        self,
        universe: List[np.ndarray],
        repcagepos: np.ndarray,
        repcagetype: List[str],
        cagetypes: Dict[str, Set[int]],
        filled_cages: Set[int],
        guests: Dict[str, Dict[str, float]],
        spot_guests: Dict[int, str],
        spot_groups: Dict[int, str],
        groups: Dict[int, Dict[int, str]],
        dopants: Dict[int, str],
        reppositions: np.ndarray,
        rotmatrices: np.ndarray,
        repcell: Cell,
    ):
        self.universe = universe
        self.repcagepos = repcagepos
        self.repcagetype = repcagetype
        self.cagetypes = cagetypes
        self.filled_cages = filled_cages
        self.guests = guests
        self.spot_guests = spot_guests
        self.spot_groups = spot_groups
        self.groups = groups
        self.dopants = dopants
        self.reppositions = reppositions
        self.rotmatrices = rotmatrices
        self.repcell = repcell

    def _add_group(self, cage, group, root):
        self.groups[root][cage] = group
        self.filled_cages.add(cage)

    @timeit
    @banner
    def execute(self) -> Stage7Output:
        """Place guest molecules."""
        logger = getLogger()

        if self.repcagepos is None:
            return Stage7Output(self.universe)

        # 1. ドーパント周辺のケージの処理
        dopants_neighbors = dopants_info(
            self.dopants, self.reppositions, self.repcagepos, self.repcell
        )

        # 2. グループの配置
        self._process_groups(dopants_neighbors, logger)

        # 3. ゲスト分子の配置
        molecules = self._process_guests(logger)

        # 4. ドーパントの処理
        self._process_dopants(logger)

        return Stage7Output(self.universe)

    def _process_groups(self, dopants_neighbors, logger):
        """グループの配置を処理"""
        if len(self.spot_groups) > 0:
            # process the -H option
            for cage, group_to in self.spot_groups.items():
                group, root = group_to.split(":")
                self._add_group(cage, group, int(root))

            logger.info("  Summary of groups:")
            groups_info(self.groups)

        # semi-guests
        for root, cages in self.groups.items():
            assert root in self.dopants
            name = self.dopants[root]
            molname = f"G{root}"
            pos = np.array(self.reppositions[root])
            rot = self.rotmatrices[root]
            self.universe.append(AtomicStructure(pos, self.repcell, name))
            del self.dopants[root]  # processed.

            for cage, group in cages.items():
                assert cage in dopants_neighbors[root]
                cage_center = self.repcagepos[cage]
                self.universe.append(
                    Group(group).arrange_atoms(
                        cage_center, pos, self.repcell, molname, origin_atom=name
                    )
                )

    def _process_guests(self, logger):
        """ゲスト分子の配置を処理"""
        molecules = defaultdict(list)

        # 特定のケージにゲストを配置
        if len(self.spot_guests) > 0:
            for cage, molec in self.spot_guests.items():
                molecules[molec].append(cage)
                self.filled_cages.add(cage)

        # ケージタイプごとのゲスト配置
        for cagetype, contents in self.guests.items():
            if cagetype not in self.cagetypes:
                logger.info(f"  Nonexistent cage type: {cagetype}")
                continue

            resident = dict()
            rooms = list(self.cagetypes[cagetype] - self.filled_cages)
            for room in rooms:
                resident[room] = None

            vacant = len(rooms)
            for molec, frac in contents.items():
                nmolec = int(frac * len(rooms) + 0.5)
                vacant -= nmolec
                assert vacant >= 0, "Too many guests."
                remain = nmolec
                movedin = []

                while remain > 0:
                    r = np.random.randint(0, len(rooms))
                    room = rooms[r]
                    if resident[room] is None:
                        resident[room] = molec
                        molecules[molec].append(room)
                        movedin.append(room)
                        remain -= 1

        # ゲスト分子の配置情報を表示
        if len(molecules):
            logger.info("  Summary of guest placements:")
            guests_info(self.cagetypes, molecules)

        # 分子ゲストの配置
        for molec, cages in molecules.items():
            guest_type, guest_options = plugin_option_parser(molec)
            logger.debug(f"  Guest type: {guest_type}")
            gmol = safe_import("molecule", guest_type).Molecule(**guest_options)

            try:
                mdoc = gmol.__doc__.splitlines()
            except BaseException:
                mdoc = []
            for line in mdoc:
                logger.info("  " + line)

            cage_center = np.array([self.repcagepos[i] for i in cages])
            cmat = [np.identity(3) for i in cages]
            self.universe.append(AtomicStructure(cage_center, self.repcell, cmat, gmol))

        return molecules

    def _process_dopants(self, logger):
        """ドーパントの処理"""
        if len(self.dopants):
            logger.info(f"  Dopants: {self.dopants}")

        # ドーパントは単原子で、1つの水分子を置換すると仮定
        atomset = defaultdict(set)
        for index, name in self.dopants.items():
            atomset[name].add(index)

        for name, indices in atomset.items():
            pos = self.reppositions[sorted(indices)]
            self.universe.append(AtomicStructure(pos, self.repcell, name=name))
