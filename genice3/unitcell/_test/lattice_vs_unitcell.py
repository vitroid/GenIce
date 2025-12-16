# 分子の向きは違えど、おおよそ同じ形であることを検証する。

import sys
import numpy as np
from logging import basicConfig, INFO, getLogger, DEBUG
import networkx as nx
from dataclasses import dataclass


@dataclass
class Product:
    oxygens: np.ndarray
    graph: nx.Graph
    volume: float
    cell: np.ndarray


def genice2_generate_ice(ice, options):

    from genice2.genice import GenIce
    from genice2.plugin import Lattice, Format, Molecule

    lattice = Lattice(ice, **options)
    formatter = Format(
        "raw",
        stage=(
            1,  # replicated molecular positions
            2,  # replicated graph
            6,  # atomic positions of water
        ),
    )
    water = Molecule("spce")
    genice2 = GenIce(lattice).generate_ice(formatter, water=water)
    cell = genice2["repcell"]
    volume = np.linalg.det(cell)
    oxygens = np.array([sites[0, :] for sites in genice2["mols"][0].positions])
    return Product(oxygens=oxygens, graph=genice2["graph"], volume=volume, cell=cell)


def GenIce3_generate_ice(ice, options):
    from genice3.plugin import UnitCell
    from genice3.plugin import Molecule
    from genice3.genice import GenIce3
    from genice3.util import validate_ice_rules

    genice3 = GenIce3()
    genice3.unitcell = UnitCell(ice, **options)

    if not validate_ice_rules(genice3.digraph):
        # ice rulesに従わない
        raise ValueError(f"Violation of the ice rules in {ice} {options=}")

    waters = genice3.water_molecules(water_model=Molecule("spce"))
    cell = genice3.cell
    volume = np.linalg.det(cell)
    oxygens = np.array([mol.sites[0] for mol in waters.values()])
    # assert the volumes are similar
    return Product(oxygens=oxygens, graph=genice3.graph, volume=volume, cell=cell)


basicConfig(level=INFO)
logger = getLogger()
for ice in sys.argv[1:]:

    if ice in ("iceMd", "iceT"):
        continue
    if ice == "xFAU":
        options = [dict(rep=4), dict(rep=8)]
    elif ice in ("one", "eleven"):
        options = [dict(layers="ccchchc"), dict(layers="ccchh")]
    elif ice == "11i":
        options = [dict(type=x) for x in range(1, 17)]
    elif ice == "bilayer":
        options = [dict(sw=0.1)]
    else:
        options = [dict()]

    for op in options:
        product2 = genice2_generate_ice(ice, op)
        product3 = GenIce3_generate_ice(ice, op)

        assert (
            abs(product2.volume - product3.volume) / product2.volume < 1e-6
        ), f"{product2.volume=}, {product3.volume=}"
        # 酸素原子位置を比較
        # 酸素の一対一対応を確認する(順序が変わっているかもしれないので。)
        import pairlist as pl

        pairs = [
            d
            for _, _, d in pl.pairs_iter(
                product2.oxygens,
                maxdist=0.1,
                cell=product2.cell,
                pos2=product3.oxygens,
                distance=True,
            )
            if d < 0.015
        ]
        logger.info(f"{min(pairs)=}, {max(pairs)=}")
        assert len(pairs) == len(
            product2.oxygens
        ), f"{len(pairs)=}, {len(product2.oxygens)=}"

        # これがなかなかに時間を食う。
        if "CRN" not in ice:
            assert nx.is_isomorphic(
                product2.graph, product3.graph
            ), f"{product2.graph.size()=}, {product3.graph.size()=}"

    # 水素位置の照合はしない。
    # bond_centers = []
    # for i, j in genice3.graph.edges():
    #     d = relative_O3[j] - relative_O3[i]
    #     d -= np.floor(d + 0.5)
    #     bond_centers.append(d / 2 + relative_O3[i])
    # bond_centers = np.array(bond_centers)

    # pairs = [
    #     d
    #     for x, y, d in pl.pairs_iter(
    #         relative_H2, maxdist=0.1, cell=cell2, pos2=bond_centers, distance=True
    #     )
    # ]
    # print(f"{min(pairs)=}, {max(pairs)=}")
    # assert len(pairs) == len(relative_H2), f"{pairs=}"

    # pairs = [
    #     d
    #     for x, y, d in pl.pairs_iter(
    #         relative_H3, maxdist=0.093, cell=cell2, pos2=bond_centers, distance=True
    #     )
    # ]
    # print(f"{min(pairs)=}, {max(pairs)=}")
    # assert len(pairs) == len(relative_H3), f"{len(pairs)=}, {len(relative_H3)=}"
