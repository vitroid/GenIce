# It is a dummy unitcell to load a cif file.
# usage: genice3 CIF[file=diamond.cif, site=C]


import genice3.unitcell
import numpy as np
from typing import Dict, List, Any, Tuple
import networkx as nx
from cif2ice import cellvectors, read_cif
from genice3.util import atoms_to_waters, shortest_distance
import re
from logging import getLogger

from genice3.cli.pool_parser import (
    parse_options_generic,
    OPTION_TYPE_STRING,
)


def parse_options(options: Dict[str, Any]) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    """
    unitcell.cifプラグインのオプションを処理

    この関数は、動的プラグインチェーン実行システムから呼び出されます。
    コマンドライン引数や設定ファイルから受け取ったオプションのうち、
    cifプラグインが処理すべきオプション（osite, file, water_model）を
    処理し、それ以外は次のプラグイン（例: moleculeプラグイン）に渡すために返します。

    Args:
        options: オプションの辞書（設定ファイルの値が初期値として含まれる可能性がある）
            - osite: 酸素位置に読みかえる原子の名前 (default "O")
            - hsite: 酸素位置に読みかえる原子の名前 (default None)
            - file: CIFファイルのパス
            - water_model: 水分子モデル名（例: "3site", "foursite"）

    Returns:
        (処理したオプション, 処理しなかったオプション) のタプル
        - 処理したオプション: gromacsプラグインが処理したオプション（osite, file, water_model）
        - 処理しなかったオプション: 次のプラグインに渡すオプション

    注意:
        - water_modelが"foursite"などのmoleculeプラグイン名の場合、
          チェーン実行システムが自動的に該当するmoleculeプラグインを呼び出します。
    """
    # オプションの型定義
    option_specs = {
        "osite": OPTION_TYPE_STRING,
        "hsite": OPTION_TYPE_STRING,
        "file": OPTION_TYPE_STRING,  # "0=foursite" または {"0": "foursite"} 形式
        "water_model": OPTION_TYPE_STRING,  # "3site", "foursite" など
    }

    # parse_options_genericを使用してオプションを処理
    # これにより、guestとspot_guestは辞書形式に変換され、
    # water_modelは文字列として処理されます
    return parse_options_generic(options, option_specs)


class UnitCell(genice3.unitcell.UnitCell):
    """
    cif単位胞を定義するクラス。
    """

    def __init__(self, **kwargs):
        logger = getLogger("CIF")
        file = kwargs.get("file")
        osite = kwargs.get("osite")
        hsite = kwargs.get("hsite")
        if file is None:
            raise ValueError("file is required")
        if osite is None:
            osite = "O"

        # download(URL, fNameIn)

        atoms, box = read_cif.read_and_process(file, make_rect_box=False)
        # pattern matching
        oatoms = np.array([a[1:] for a in atoms if re.match(osite, a[0])])
        logger.info(f"{osite=} {oatoms.shape=} {atoms}")
        cell = cellvectors(*box) / 10  # nm
        shortest_OO = shortest_distance(oatoms, cell)
        cell *= 0.276 / shortest_OO

        if hsite is None:
            # 水素位置は指定されていないので、genice3にまかせる。
            super().__init__(
                cell=cell,  # nm
                lattice_sites=oatoms,
                coord="relative",
            )
        else:
            hatoms = np.array([a[1:] for a in atoms if re.match(hsite, a[0])])
            waters, pairs, oo_pairs = atoms_to_waters(
                oatoms, hatoms, cell, partial_order=True
            )
            super().__init__(
                cell=cell,
                lattice_sites=waters,
                graph=nx.Graph(oo_pairs),
                fixed=nx.DiGraph(pairs),
                coord="relative",
            )
