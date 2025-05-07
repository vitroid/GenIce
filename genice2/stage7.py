from dataclasses import dataclass
from typing import List, Dict, Set
from logging import getLogger
from collections import defaultdict

import numpy as np

from genice2.decorators import banner, timeit
from genice2.molecules import Molecule, arrange, one
from genice2.plugin import Group, safe_import
from genice2.valueparser import plugin_option_parser
from genice2.cell import Cell, rel_wrap


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
                # logger.info((i,cagepos[i]))

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

    @timeit
    @banner
    def execute(self) -> Stage7Output:
        """Place guest molecules."""
        # ゲスト分子の配置処理

        logger = getLogger()

        if self.repcagepos is not None:
            # the cages around the dopants.
            dopants_neighbors = dopants_info(
                self.dopants, self.reppositions, self.repcagepos, self.repcell
            )

            # put the (one-off) groups
            if len(self.spot_groups) > 0:
                # process the -H option
                for cage, group_to in self.spot_groups.items():
                    group, root = group_to.split(":")
                    self.add_group(cage, group, int(root))

            molecules = defaultdict(list)

            if len(self.spot_guests) > 0:
                # process the -G option
                for cage, molec in self.spot_guests.items():
                    molecules[molec].append(cage)
                    self.filled_cages.add(cage)

            # process the -g option
            for cagetype, contents in self.guests.items():
                if cagetype not in self.cagetypes:
                    logger.info(f"Nonexistent cage type: {cagetype}")
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

            # Now ge got the address book of the molecules.
            if len(molecules):
                logger.info("  Summary of guest placements:")
                guests_info(self.cagetypes, molecules)

            if len(self.spot_groups) > 0:
                logger.info("  Summary of groups:")
                groups_info(self.groups)

            # semi-guests
            for root, cages in self.groups.items():
                assert root in self.dopants
                name = self.dopants[root]
                molname = f"G{root}"
                pos = self.reppositions[root]
                rot = self.rotmatrices[root]
                self.universe.append(monatom(pos, self.repcell, name))
                del self.dopants[root]  # processed.
                logger.debug((root, cages, name, molname, pos, rot))

                for cage, group in cages.items():
                    # assert group in self.groups_placer
                    assert cage in dopants_neighbors[root]
                    cage_center = self.repcagepos[cage]
                    self.universe.append(
                        Group(group).arrange_atoms(
                            cage_center, pos, self.repcell, molname, origin_atom=name
                        )
                    )

            # molecular guests
            for molec, cages in molecules.items():
                guest_type, guest_options = plugin_option_parser(molec)
                logger.debug(f"Guest type: {guest_type}")
                gmol = safe_import("molecule", guest_type).Molecule(**guest_options)

                try:
                    mdoc = gmol.__doc__.splitlines()
                except BaseException:
                    mdoc = []
                for line in mdoc:
                    logger.info("  " + line)
                cage_center = [self.repcagepos[i] for i in cages]
                cmat = [np.identity(3) for i in cages]
                self.universe.append(arrange(cage_center, self.repcell, cmat, gmol))

        logger.info(f"  Dopants: {self.dopants}")

        # Assume the dopant is monatomic and replaces one water molecule
        atomset = defaultdict(set)
        for label, name in self.dopants.items():
            atomset[name].add(label)

        for name, labels in atomset.items():
            pos = [self.reppositions[i] for i in sorted(labels)]
            rot = [self.rotmatrices[i] for i in sorted(labels)]
            oneatom = one.Molecule(label=name)
            self.universe.append(arrange(pos, self.repcell, rot, oneatom))

        return Stage7Output(self.universe)
