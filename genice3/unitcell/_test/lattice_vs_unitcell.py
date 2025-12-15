# 分子の向きは違えど、おおよそ同じ形であることを検証する。

import sys
import numpy as np
from logging import basicConfig, INFO, getLogger, DEBUG
import networkx as nx

basicConfig(level=DEBUG)
logger = getLogger()
for ice in sys.argv[1:]:

    from genice2.genice import GenIce
    from genice2.plugin import Lattice, Format, Molecule

    if ice in ("iceMd", "iceT"):
        continue
    if ice == "xFAU":
        lattice = Lattice(ice, rep=4)
    elif ice == "one":
        lattice = Lattice(ice, layers="ccchchc")
    else:
        lattice = Lattice(ice)
    formatter = Format(
        "raw",
        stage=(
            1,
            2,
            6,
        ),
    )
    water = Molecule("spce")
    genice2 = GenIce(lattice).generate_ice(formatter, water=water)
    cell2 = genice2["repcell"]
    volume2 = np.linalg.det(cell2)
    cell2i = np.linalg.inv(cell2)
    O2 = np.array([sites[0, :] for sites in genice2["mols"][0].positions])

    import numpy as np

    H2 = np.array(
        [sites[1, :] for sites in genice2["mols"][0].positions]
        + [sites[2, :] for sites in genice2["mols"][0].positions]
    )
    relative_O2 = O2 @ cell2i
    relative_H2 = H2 @ cell2i
    assert len(relative_O2) * 2 == len(relative_H2), f"{relative_O2=}, {relative_H2=}"

    from genice3.plugin import UnitCell
    from genice3.plugin import Molecule
    from genice3.genice import GenIce3

    genice3 = GenIce3()
    if ice == "xFAU":
        genice3.unitcell = UnitCell(ice, rep=4)
    elif ice == "one":
        genice3.unitcell = UnitCell(ice, layers="ccchchc")
    else:
        genice3.unitcell = UnitCell(ice)
    waters = genice3.water_molecules(water_model=Molecule("spce"))
    cell3 = genice3.unitcell.cell
    volume3 = np.linalg.det(cell3)
    # assert the volumes are similar
    assert abs(volume2 - volume3) / volume2 < 1e-6, f"{volume2=}, {volume3=}"
    cell3i = np.linalg.inv(cell3)
    O3 = [mol.sites[0] for mol in waters.values()]
    H3 = np.vstack(
        [
            [mol.sites[1] for mol in waters.values()],
            [mol.sites[2] for mol in waters.values()],
        ]
    )
    relative_O3 = O3 @ cell3i
    relative_H3 = H3 @ cell3i
    assert len(relative_O3) * 2 == len(relative_H3), f"{relative_O3=}, {relative_H3=}"
    # 酸素原子位置を比較
    # 酸素の一対一対応を確認する(順序が変わっているかもしれないので。)
    import pairlist as pl

    pairs = [
        (x, y, d)
        for x, y, d in pl.pairs_iter(
            relative_O2, maxdist=0.1, cell=cell2, pos2=relative_O3, distance=True
        )
        if d < 0.015
    ]
    assert len(pairs) == len(relative_O2)

    # print(genice2)
    assert nx.is_isomorphic(genice2["graph"], genice3.graph)

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
