"""
GenIce3をDagsterのassetに変換した実装例。

設計方針:
1. GenIce3クラスは解体し、すべてを関数とassetに変換
2. 設定パラメータは@dataclassで管理
3. 計算結果はDagsterのassetとして実装
4. 依存関係はasset依存関係として宣言

この設計により、関数型プログラミングの原則に従い、
内部状態を持たない純粋な関数として実装できる。
"""

from dataclasses import dataclass
from typing import Dict, Tuple, Optional
import numpy as np
import networkx as nx
from dagster import asset, Definitions
from genice2.unitcell import UnitCell
from genice2.molecules import Molecule
from genice2 import ConfigurationError
from genice2.stage1 import replicate_positions
from genice2.stage2 import grandcell_wrap
from genice2.cell import cellvectors, cellshape
from genice2.stage5 import assume_tetrahedral_vectors
import genice_core
from genice2.genice3 import (
    orientations as compute_orientations,
    replicate_graph,
    replicate_fixed_edges,
)


@dataclass
class GenIce3Config:
    """
    GenIce3の設定パラメータを管理するデータクラス。

    クラス変数のような内部状態ではなく、明示的な設定オブジェクトとして管理する。
    """

    depol_loop: int = 1000
    replication_matrix: np.ndarray = None

    def __post_init__(self):
        if self.replication_matrix is None:
            self.replication_matrix = np.eye(3, dtype=int)


# ============================================================================
# Asset定義: 各計算ステップを独立したassetとして実装
# ============================================================================


# 設定と単位胞をassetとして定義（Configを使わず、シンプルにデフォルト値で作成）
@asset
def genice3_config() -> GenIce3Config:
    """
    GenIce3の設定を生成するasset。
    デフォルト値を使用します。
    """
    return GenIce3Config(
        depol_loop=1000,
        replication_matrix=np.eye(3, dtype=int),
    )


@asset
def unitcell_asset() -> UnitCell:
    """
    単位胞を生成するasset。
    デフォルト値（A15）を使用します。
    """
    from genice2.unitcell import A15

    return A15(
        shift=(0.0, 0.0, 0.0),
        assess_cages=False,
        anions={},
        cations={},
    )


@asset
def cell(
    genice3_config: GenIce3Config,
    unitcell_asset: UnitCell,
) -> np.ndarray:
    """
    単位胞と複製行列から拡大単位胞のセルベクトルを計算する。

    Args:
        genice3_config: GenIce3の設定
        unitcell_asset: 単位胞オブジェクト

    Returns:
        拡大単位胞のセルベクトル
    """
    return unitcell_asset.cell @ genice3_config.replication_matrix


@asset
def replica_vectors(genice3_config: GenIce3Config) -> np.ndarray:
    """
    レプリカベクトルの計算。

    Args:
        config: GenIce3の設定

    Returns:
        レプリカベクトルの配列
    """
    i, j, k = np.array(genice3_config.replication_matrix)
    corners = np.array(
        [a * i + b * j + c * k for a in (0, 1) for b in (0, 1) for c in (0, 1)]
    )

    mins = np.min(corners, axis=0)
    maxs = np.max(corners, axis=0)

    det = abs(np.linalg.det(genice3_config.replication_matrix))
    det = np.floor(det + 0.5).astype(int)
    invdet = np.floor(
        np.linalg.inv(genice3_config.replication_matrix) * det + 0.5
    ).astype(int)

    vecs = set()
    for a in range(mins[0], maxs[0] + 1):
        for b in range(mins[1], maxs[1] + 1):
            for c in range(mins[2], maxs[2] + 1):
                abc = np.array([a, b, c])
                rep = grandcell_wrap(
                    abc, genice3_config.replication_matrix, invdet, det
                ).astype(int)
                if tuple(rep) not in vecs:
                    vecs.add(tuple(rep))

    vecs = np.array(list(vecs))
    vol = abs(np.linalg.det(genice3_config.replication_matrix))
    assert np.allclose(vol, len(vecs)), (vol, vecs)

    return vecs


@asset
def replica_vector_labels(
    replica_vectors: np.ndarray,
) -> Dict[Tuple[int, int, int], int]:
    """
    レプリカベクトル座標をラベルにマッピングする辞書を生成。

    Args:
        replica_vectors: レプリカベクトルの配列

    Returns:
        レプリカベクトル座標をキー、ラベルを値とする辞書
    """
    # numpyのint64をPythonのintに変換
    return {tuple(int(x) for x in xyz): int(i) for i, xyz in enumerate(replica_vectors)}


@asset
def lattice_sites(
    unitcell_asset: UnitCell,
    replica_vectors: np.ndarray,
    genice3_config: GenIce3Config,
) -> np.ndarray:
    """
    格子点の位置を複製する。

    Args:
        unitcell_asset: 単位胞オブジェクト
        replica_vectors: レプリカベクトル
        genice3_config: GenIce3の設定

    Returns:
        複製された格子点の位置
    """
    return replicate_positions(
        unitcell_asset.waters, replica_vectors, genice3_config.replication_matrix
    )


@asset
def graph(
    unitcell_asset: UnitCell,
    replica_vectors: np.ndarray,
    replica_vector_labels: Dict[Tuple[int, int, int], int],
    genice3_config: GenIce3Config,
) -> nx.Graph:
    """
    グラフを複製する。

    Args:
        unitcell_asset: 単位胞オブジェクト
        replica_vectors: レプリカベクトル
        replica_vector_labels: レプリカベクトルラベル
        genice3_config: GenIce3の設定

    Returns:
        複製されたグラフ
    """
    return replicate_graph(
        unitcell_asset.graph,
        unitcell_asset.waters,
        replica_vectors,
        replica_vector_labels,
        genice3_config.replication_matrix,
    )


@asset
def fixed_edges(
    graph: nx.Graph,
    unitcell_asset: UnitCell,
) -> nx.DiGraph:
    """
    固定エッジを複製する。

    Args:
        graph: 複製されたグラフ
        unitcell_asset: 単位胞オブジェクト

    Returns:
        複製された固定エッジ
    """
    return replicate_fixed_edges(
        graph, unitcell_asset.fixed, len(unitcell_asset.waters)
    )


@asset
def digraph(
    graph: nx.Graph,
    lattice_sites: np.ndarray,
    fixed_edges: nx.DiGraph,
    genice3_config: GenIce3Config,
) -> nx.DiGraph:
    """
    有向グラフを生成する。

    Args:
        graph: 複製されたグラフ
        lattice_sites: 格子点の位置
        fixed_edges: 固定エッジ
        genice3_config: GenIce3の設定

    Returns:
        有向グラフ
    """
    dg = genice_core.ice_graph(
        graph,
        vertexPositions=lattice_sites,
        isPeriodicBoundary=True,
        dipoleOptimizationCycles=genice3_config.depol_loop,
        fixedEdges=fixed_edges,
    )
    if not dg:
        raise ConfigurationError("Failed to generate a directed graph.")
    return dg


@asset
def anions(
    unitcell_asset: UnitCell,
    replica_vectors: np.ndarray,
) -> Dict[int, str]:
    """
    アニオン位置を複製する。

    Args:
        unitcell_asset: 単位胞オブジェクト
        replica_vectors: レプリカベクトル

    Returns:
        アニオン位置の辞書
    """
    anion_dict = dict()
    Z = len(unitcell_asset.waters)
    for label, ion_name in unitcell_asset.anions.items():
        for i in range(len(replica_vectors)):
            site = i * Z + label
            anion_dict[site] = ion_name
    return anion_dict


@asset
def cations(
    unitcell_asset: UnitCell,
    replica_vectors: np.ndarray,
) -> Dict[int, str]:
    """
    カチオン位置を複製する。

    Args:
        unitcell_asset: 単位胞オブジェクト
        replica_vectors: レプリカベクトル

    Returns:
        カチオン位置の辞書
    """
    cation_dict = dict()
    Z = len(unitcell_asset.waters)
    for label, ion_name in unitcell_asset.cations.items():
        for i in range(len(replica_vectors)):
            site = i * Z + label
            cation_dict[site] = ion_name
    return cation_dict


@asset
def site_occupants(
    lattice_sites: np.ndarray,
    anions: Dict[int, str],
    cations: Dict[int, str],
) -> list:
    """
    格子点の占有者を決定する。

    Args:
        lattice_sites: 格子点の位置
        anions: アニオン位置
        cations: カチオン位置

    Returns:
        各格子点の占有者のリスト
    """
    occupants = ["water"] * len(lattice_sites)
    for site, ion_name in anions.items():
        occupants[site] = ion_name
    for site, ion_name in cations.items():
        occupants[site] = ion_name
    return occupants


@asset
def molecule_orientations(
    lattice_sites: np.ndarray,
    digraph: nx.DiGraph,
    cell: np.ndarray,
    anions: Dict[int, str],
    cations: Dict[int, str],
) -> np.ndarray:
    """
    水分子の向きを計算する。

    Args:
        lattice_sites: 格子点の位置
        digraph: 有向グラフ
        cell: セルベクトル
        anions: アニオン位置
        cations: カチオン位置

    Returns:
        各格子点の回転行列
    """
    return compute_orientations(
        lattice_sites,
        digraph,
        cell,
        anions | cations,
    )


@asset
def cages(
    unitcell_asset: UnitCell,
    replica_vectors: np.ndarray,
    genice3_config: GenIce3Config,
) -> Tuple[np.ndarray, list]:
    """
    ケージ位置を複製する。

    Args:
        unitcell_asset: 単位胞オブジェクト
        replica_vectors: レプリカベクトル
        genice3_config: GenIce3の設定

    Returns:
        ケージ位置とタイプのタプル
    """
    if unitcell_asset.cages is None:
        raise ValueError("Unit cell does not have cages defined.")
    repcagepos = replicate_positions(
        unitcell_asset.cages[0], replica_vectors, genice3_config.replication_matrix
    )
    num_cages_in_unitcell = len(unitcell_asset.cages[0])
    repcagetype = [
        unitcell_asset.cages[1][i % num_cages_in_unitcell]
        for i in range(len(repcagepos))
    ]
    return repcagepos, repcagetype


# ============================================================================
# 分子生成のasset
# ============================================================================


@asset
def waters(
    lattice_sites: np.ndarray,
    molecule_orientations: np.ndarray,
    cell: np.ndarray,
    site_occupants: list,
) -> Dict[int, Molecule]:
    """
    水分子を生成するasset。

    Args:
        lattice_sites: 格子点の位置
        molecule_orientations: 回転行列
        cell: セルベクトル
        site_occupants: 格子点の占有者

    Returns:
        水分子の辞書
    """
    from genice2.genice3 import FourSiteWater

    water_model = FourSiteWater()
    mols = {}
    for site in range(len(lattice_sites)):
        if "water" == site_occupants[site]:
            rel_position = lattice_sites[site]
            orientation = molecule_orientations[site]

            sites = water_model.sites @ orientation + rel_position @ cell
            mols[site] = Molecule(
                name=water_model.name,
                sites=sites,
                labels=water_model.labels,
                is_water=True,
            )
    return mols


@asset
def ions(
    lattice_sites: np.ndarray,
    cell: np.ndarray,
    anions: Dict[int, str],
    cations: Dict[int, str],
) -> Dict[int, Molecule]:
    """
    置換イオンを生成するasset。

    Args:
        lattice_sites: 格子点の位置
        cell: セルベクトル
        anions: アニオン位置
        cations: カチオン位置

    Returns:
        イオンの辞書
    """
    ion_dict = {}
    for label, name in anions.items():
        ion_dict[label] = Molecule(
            name=name,
            sites=[lattice_sites[label] @ cell],
            labels=[name],
            is_water=False,
        )
    for label, name in cations.items():
        ion_dict[label] = Molecule(
            name=name,
            sites=[lattice_sites[label] @ cell],
            labels=[name],
            is_water=False,
        )
    return ion_dict


@asset
def guests(
    cages: Tuple[np.ndarray, list],
    cell: np.ndarray,
) -> list:
    """
    ゲスト分子を生成するasset。

    デフォルトでは空の辞書を使用するため、結果として空のリストを返します。
    CLIではパラメータに基づいて動的に生成されます。

    Args:
        cages: ケージ位置とタイプ
        cell: セルベクトル

    Returns:
        ゲスト分子のリスト
    """
    # デフォルトでは空の辞書を使用（ゲスト分子なし）
    # これにより、ionsと同様に実際の処理ロジックを通る
    return guest_molecules(cages, cell, {})


@asset
def gro_output(
    cell: np.ndarray,
    waters: Dict[int, Molecule],
    ions: Dict[int, Molecule],
    guests: list,
) -> str:
    """
    Gromacs形式のgroファイルの内容を生成するasset。

    Args:
        cell: セルベクトル
        waters: 水分子の辞書
        ions: イオンの辞書
        guests: ゲスト分子のリスト

    Returns:
        Gromacs形式の文字列
    """
    return to_gro(cell, waters, guests, ions)


# ============================================================================
# ユーティリティ関数: メソッドを関数に変換（CLI用）
# ============================================================================


def water_molecules(
    lattice_sites: np.ndarray,
    orientations: np.ndarray,
    cell: np.ndarray,
    site_occupants: list,
    water_model: Molecule,
) -> Dict[int, Molecule]:
    """
    水分子を生成する。

    Args:
        lattice_sites: 格子点の位置
        orientations: 回転行列
        cell: セルベクトル
        site_occupants: 格子点の占有者
        water_model: 水分子モデル

    Returns:
        水分子の辞書
    """
    mols = {}
    for site in range(len(lattice_sites)):
        if "water" == site_occupants[site]:
            rel_position = lattice_sites[site]
            orientation = orientations[site]

            sites = water_model.sites @ orientation + rel_position @ cell
            mols[site] = Molecule(
                name=water_model.name,
                sites=sites,
                labels=water_model.labels,
                is_water=True,
            )
    return mols


def guest_molecules(
    cages: Tuple[np.ndarray, list],
    cell: np.ndarray,
    guests: Dict[str, list],  # GuestSpecのリスト
    seed: Optional[int] = None,
) -> list:
    """
    ゲスト分子を生成する。

    Args:
        cages: ケージ位置とタイプ
        cell: セルベクトル
        guests: ゲスト仕様の辞書
        seed: 乱数シード（オプション）

    Returns:
        ゲスト分子のリスト
    """
    if seed is not None:
        np.random.seed(seed)

    randoms = np.random.random(len(cages[0]))

    mols = []
    for pos, label, probability in zip(cages[0], cages[1], randoms):
        accum = 0.0
        if label in guests:
            for guest_spec in guests[label]:
                molecule = guest_spec.molecule
                occupancy = guest_spec.occupancy
                accum += occupancy
                if accum > probability:
                    mols.append(
                        Molecule(
                            name=molecule.name,
                            sites=molecule.sites + pos @ cell,
                            labels=molecule.labels,
                            is_water=molecule.is_water,
                        )
                    )
                    break
    return mols


def substitutional_ions(
    lattice_sites: np.ndarray,
    cell: np.ndarray,
    anions: Dict[int, str],
    cations: Dict[int, str],
) -> Dict[int, Molecule]:
    """
    置換イオンを生成する。

    Args:
        lattice_sites: 格子点の位置
        cell: セルベクトル
        anions: アニオン位置
        cations: カチオン位置

    Returns:
        イオンの辞書
    """
    ions = {}
    for label, name in anions.items():
        ions[label] = Molecule(
            name=name,
            sites=[lattice_sites[label] @ cell],
            labels=[name],
            is_water=False,
        )
    for label, name in cations.items():
        ions[label] = Molecule(
            name=name,
            sites=[lattice_sites[label] @ cell],
            labels=[name],
            is_water=False,
        )
    return ions


# ============================================================================
# Definitions: Dagsterの定義オブジェクト
# ============================================================================

defs = Definitions(
    assets=[
        genice3_config,
        unitcell_asset,
        cell,
        replica_vectors,
        replica_vector_labels,
        lattice_sites,
        graph,
        fixed_edges,
        digraph,
        anions,
        cations,
        site_occupants,
        molecule_orientations,
        cages,
        waters,
        ions,
        guests,
        gro_output,
    ]
)


# ============================================================================
# CLI: コマンドラインインターフェース
# ============================================================================

import click
from logging import getLogger, DEBUG, INFO, basicConfig
from genice2.genice3 import guest_parser, ion_parser, FourSiteWater
from genice2.unitcell import A15, Ice1h
from genice2.cell import cellvectors, cellshape


def create_dynamic_assets(
    depol_loop: int = 1000,
    replication_matrix: Optional[np.ndarray] = None,
    unitcell_type: str = "A15",
    shift: Tuple[float, float, float] = (0.0, 0.0, 0.0),
    assess_cages: bool = False,
    anions: Dict[int, str] = None,
    cations: Dict[int, str] = None,
    density: Optional[float] = None,
    guests: Dict[str, list] = None,  # GuestSpecのリスト
):
    """
    パラメータに基づいて動的にasset定義を作成する。
    """
    if anions is None:
        anions = {}
    if cations is None:
        cations = {}
    if guests is None:
        guests = {}
    if replication_matrix is None:
        replication_matrix = np.eye(3, dtype=int)

    # 動的にassetを定義
    @asset
    def dynamic_genice3_config() -> GenIce3Config:
        return GenIce3Config(
            depol_loop=depol_loop,
            replication_matrix=replication_matrix,
        )

    @asset
    def dynamic_unitcell_asset() -> UnitCell:
        kwargs = {
            "shift": shift,
            "assess_cages": assess_cages,
            "anions": anions,
            "cations": cations,
        }
        if density is not None:
            kwargs["density"] = density

        if unitcell_type == "A15":
            return A15(**kwargs)
        elif unitcell_type == "Ice1h":
            return Ice1h(**kwargs)
        else:
            raise ValueError(f"Unknown unitcell type: {unitcell_type}")

    # 既存のassetを再利用（依存関係を変更）
    @asset
    def dynamic_cell(
        dynamic_genice3_config: GenIce3Config,
        dynamic_unitcell_asset: UnitCell,
    ) -> np.ndarray:
        return dynamic_unitcell_asset.cell @ dynamic_genice3_config.replication_matrix

    @asset
    def dynamic_replica_vectors(dynamic_genice3_config: GenIce3Config) -> np.ndarray:
        i, j, k = np.array(dynamic_genice3_config.replication_matrix)
        corners = np.array(
            [a * i + b * j + c * k for a in (0, 1) for b in (0, 1) for c in (0, 1)]
        )
        mins = np.min(corners, axis=0)
        maxs = np.max(corners, axis=0)
        det = abs(np.linalg.det(dynamic_genice3_config.replication_matrix))
        det = np.floor(det + 0.5).astype(int)
        invdet = np.floor(
            np.linalg.inv(dynamic_genice3_config.replication_matrix) * det + 0.5
        ).astype(int)
        vecs = set()
        for a in range(mins[0], maxs[0] + 1):
            for b in range(mins[1], maxs[1] + 1):
                for c in range(mins[2], maxs[2] + 1):
                    abc = np.array([a, b, c])
                    rep = grandcell_wrap(
                        abc, dynamic_genice3_config.replication_matrix, invdet, det
                    ).astype(int)
                    if tuple(rep) not in vecs:
                        vecs.add(tuple(rep))
        vecs = np.array(list(vecs))
        vol = abs(np.linalg.det(dynamic_genice3_config.replication_matrix))
        assert np.allclose(vol, len(vecs)), (vol, vecs)
        return vecs

    @asset
    def dynamic_replica_vector_labels(
        dynamic_replica_vectors: np.ndarray,
    ) -> Dict[Tuple[int, int, int], int]:
        # numpyのint64をPythonのintに変換
        return {
            tuple(int(x) for x in xyz): int(i)
            for i, xyz in enumerate(dynamic_replica_vectors)
        }

    @asset
    def dynamic_lattice_sites(
        dynamic_unitcell_asset: UnitCell,
        dynamic_replica_vectors: np.ndarray,
        dynamic_genice3_config: GenIce3Config,
    ) -> np.ndarray:
        return replicate_positions(
            dynamic_unitcell_asset.waters,
            dynamic_replica_vectors,
            dynamic_genice3_config.replication_matrix,
        )

    @asset
    def dynamic_graph(
        dynamic_unitcell_asset: UnitCell,
        dynamic_replica_vectors: np.ndarray,
        dynamic_replica_vector_labels: Dict[Tuple[int, int, int], int],
        dynamic_genice3_config: GenIce3Config,
    ) -> nx.Graph:
        return replicate_graph(
            dynamic_unitcell_asset.graph,
            dynamic_unitcell_asset.waters,
            dynamic_replica_vectors,
            dynamic_replica_vector_labels,
            dynamic_genice3_config.replication_matrix,
        )

    @asset
    def dynamic_fixed_edges(
        dynamic_graph: nx.Graph,
        dynamic_unitcell_asset: UnitCell,
    ) -> nx.DiGraph:
        return replicate_fixed_edges(
            dynamic_graph,
            dynamic_unitcell_asset.fixed,
            len(dynamic_unitcell_asset.waters),
        )

    @asset
    def dynamic_digraph(
        dynamic_graph: nx.Graph,
        dynamic_lattice_sites: np.ndarray,
        dynamic_fixed_edges: nx.DiGraph,
        dynamic_genice3_config: GenIce3Config,
    ) -> nx.DiGraph:
        dg = genice_core.ice_graph(
            dynamic_graph,
            vertexPositions=dynamic_lattice_sites,
            isPeriodicBoundary=True,
            dipoleOptimizationCycles=dynamic_genice3_config.depol_loop,
            fixedEdges=dynamic_fixed_edges,
        )
        if not dg:
            raise ConfigurationError("Failed to generate a directed graph.")
        return dg

    @asset
    def dynamic_anions(
        dynamic_unitcell_asset: UnitCell,
        dynamic_replica_vectors: np.ndarray,
    ) -> Dict[int, str]:
        anion_dict = {}
        Z = len(dynamic_unitcell_asset.waters)
        for label, ion_name in dynamic_unitcell_asset.anions.items():
            for i in range(len(dynamic_replica_vectors)):
                site = i * Z + label
                anion_dict[site] = ion_name
        return anion_dict

    @asset
    def dynamic_cations(
        dynamic_unitcell_asset: UnitCell,
        dynamic_replica_vectors: np.ndarray,
    ) -> Dict[int, str]:
        cation_dict = {}
        Z = len(dynamic_unitcell_asset.waters)
        for label, ion_name in dynamic_unitcell_asset.cations.items():
            for i in range(len(dynamic_replica_vectors)):
                site = i * Z + label
                cation_dict[site] = ion_name
        return cation_dict

    @asset
    def dynamic_site_occupants(
        dynamic_lattice_sites: np.ndarray,
        dynamic_anions: Dict[int, str],
        dynamic_cations: Dict[int, str],
    ) -> list:
        occupants = ["water"] * len(dynamic_lattice_sites)
        for site, ion_name in dynamic_anions.items():
            occupants[site] = ion_name
        for site, ion_name in dynamic_cations.items():
            occupants[site] = ion_name
        return occupants

    @asset
    def dynamic_molecule_orientations(
        dynamic_lattice_sites: np.ndarray,
        dynamic_digraph: nx.DiGraph,
        dynamic_cell: np.ndarray,
        dynamic_anions: Dict[int, str],
        dynamic_cations: Dict[int, str],
    ) -> np.ndarray:
        return compute_orientations(
            dynamic_lattice_sites,
            dynamic_digraph,
            dynamic_cell,
            dynamic_anions | dynamic_cations,
        )

    @asset
    def dynamic_cages(
        dynamic_unitcell_asset: UnitCell,
        dynamic_replica_vectors: np.ndarray,
        dynamic_genice3_config: GenIce3Config,
    ) -> Tuple[np.ndarray, list]:
        if dynamic_unitcell_asset.cages is None:
            raise ValueError("Unit cell does not have cages defined.")
        repcagepos = replicate_positions(
            dynamic_unitcell_asset.cages[0],
            dynamic_replica_vectors,
            dynamic_genice3_config.replication_matrix,
        )
        num_cages_in_unitcell = len(dynamic_unitcell_asset.cages[0])
        repcagetype = [
            dynamic_unitcell_asset.cages[1][i % num_cages_in_unitcell]
            for i in range(len(repcagepos))
        ]
        return repcagepos, repcagetype

    @asset
    def dynamic_guests(
        dynamic_cages: Tuple[np.ndarray, list],
        dynamic_cell: np.ndarray,
    ) -> list:
        """
        ゲスト分子を生成する動的asset。
        """
        # guest_molecules関数を使用
        return guest_molecules(dynamic_cages, dynamic_cell, guests)

    @asset
    def dynamic_waters(
        dynamic_lattice_sites: np.ndarray,
        dynamic_molecule_orientations: np.ndarray,
        dynamic_cell: np.ndarray,
        dynamic_site_occupants: list,
    ) -> Dict[int, Molecule]:
        """
        水分子を生成する動的asset。
        """
        from genice2.genice3 import FourSiteWater

        return water_molecules(
            dynamic_lattice_sites,
            dynamic_molecule_orientations,
            dynamic_cell,
            dynamic_site_occupants,
            FourSiteWater(),
        )

    @asset
    def dynamic_ions(
        dynamic_lattice_sites: np.ndarray,
        dynamic_cell: np.ndarray,
        dynamic_anions: Dict[int, str],
        dynamic_cations: Dict[int, str],
    ) -> Dict[int, Molecule]:
        """
        置換イオンを生成する動的asset。
        """
        return substitutional_ions(
            dynamic_lattice_sites,
            dynamic_cell,
            dynamic_anions,
            dynamic_cations,
        )

    @asset
    def dynamic_gro_output(
        dynamic_cell: np.ndarray,
        dynamic_waters: Dict[int, Molecule],
        dynamic_ions: Dict[int, Molecule],
        dynamic_guests: list,
    ) -> str:
        """
        Gromacs形式のgroファイルの内容を生成する動的asset。
        """
        return to_gro(dynamic_cell, dynamic_waters, dynamic_guests, dynamic_ions)

    return {
        "genice3_config": dynamic_genice3_config,
        "unitcell_asset": dynamic_unitcell_asset,
        "cell": dynamic_cell,
        "replica_vectors": dynamic_replica_vectors,
        "replica_vector_labels": dynamic_replica_vector_labels,
        "lattice_sites": dynamic_lattice_sites,
        "graph": dynamic_graph,
        "fixed_edges": dynamic_fixed_edges,
        "digraph": dynamic_digraph,
        "anions": dynamic_anions,
        "cations": dynamic_cations,
        "site_occupants": dynamic_site_occupants,
        "molecule_orientations": dynamic_molecule_orientations,
        "cages": dynamic_cages,
        "guests": dynamic_guests,
        "waters": dynamic_waters,
        "ions": dynamic_ions,
        "gro_output": dynamic_gro_output,
    }


@click.command()
@click.help_option()
@click.option("--debug", "-D", is_flag=True, help="Enable debug mode")
@click.option(
    "--shift",
    "-S",
    type=click.Tuple([float, float, float]),
    default=(0.0, 0.0, 0.0),
    help="Shift the unit cell",
)
@click.option("--depol_loop", type=int, default=1000, help="Depolarization loop")
@click.option(
    "--replication_matrix",
    type=click.Tuple([int, int, int, int, int, int, int, int, int]),
    default=None,
    help="Replication matrix (9 elements)",
)
@click.option(
    "--replication_factors",
    "--rep",
    type=click.Tuple([int, int, int]),
    default=(1, 1, 1),
    help="Replication factors (3 elements)",
)
@click.option("--assess_cages", "-A", is_flag=True, help="Assess cages")
@click.option(
    "--guest", "-g", multiple=True, help="Guest descriptors, e.g. A12=me, A14=co2, etc."
)
@click.option(
    "--anion",
    "-a",
    multiple=True,
    help="Specify a monatomic anion that replaces a water molecule., e.g. -a 1=Cl, -a 35=Br, etc.",
)
@click.option(
    "--cation",
    "-c",
    multiple=True,
    help="Specify a monatomic cation that replaces a water molecule., e.g. -c 1=Na, -c 35=K, etc.",
)
@click.option(
    "--density", "--dens", type=float, default=None, help="Density of the ice"
)
@click.option(
    "--unitcell_type",
    type=click.Choice(["A15", "Ice1h"]),
    default="A15",
    help="Unit cell type",
)
@click.option(
    "--output",
    "-o",
    type=str,
    default="genice3_dagster.gro",
    help="Output file name",
)
def main(
    debug,
    shift,
    depol_loop,
    replication_matrix,
    replication_factors,
    assess_cages,
    guest,
    anion,
    cation,
    density,
    unitcell_type,
    output,
):
    """
    GenIce3 Dagster版のCLI。
    Dagsterのassetを実行して、結果をgroファイルに出力します。
    """
    basicConfig(level=DEBUG if debug else INFO)
    logger = getLogger()

    # パラメータの処理
    if replication_matrix is None:
        replication_matrix = np.diag(replication_factors)
    else:
        replication_matrix = np.array(replication_matrix).reshape(3, 3)

    guest_info = guest_parser(guest) if guest else {}
    anion_info = ion_parser(anion) if anion else {}
    cation_info = ion_parser(cation) if cation else {}

    logger.info("Creating dynamic assets...")
    # 動的にassetを作成
    dynamic_assets = create_dynamic_assets(
        depol_loop=depol_loop,
        replication_matrix=replication_matrix,
        unitcell_type=unitcell_type,
        shift=shift,
        assess_cages=assess_cages,
        anions=anion_info,
        cations=cation_info,
        density=density,
        guests=guest_info,
    )

    logger.info("Materializing assets...")
    # assetを実行（依存関係に基づいて自動的に必要なassetが実行される）
    from dagster import materialize

    # すべてのassetを渡すことで、依存関係に基づいて自動的に実行される
    # 最終的なassetだけを指定しても、依存関係にあるすべてのassetが実行される
    result = materialize(
        list(dynamic_assets.values()),
    )

    # 結果を取得（gro_output assetから直接取得）
    logger.info(f"Writing output to {output}...")
    gro_content = result.output_for_node("dynamic_gro_output")

    with open(output, "w") as f:
        f.write(gro_content)

    logger.info(f"Successfully generated {output}")


def to_gro(cellmat: np.ndarray, waters: dict, guests: list, ions: dict) -> str:
    """
    Gromacs形式で出力する。

    Args:
        cellmat: セル行列
        waters: 水分子の辞書
        guests: ゲスト分子のリスト
        ions: イオンの辞書

    Returns:
        Gromacs形式の文字列
    """
    logger = getLogger()

    if not (cellmat[0, 1] == 0 and cellmat[0, 2] == 0 and cellmat[1, 2] == 0):
        logger.info(
            "  The specified reshaping matrix does not obey the requirements for Gromacs' unit cell convention."
        )
        a = np.linalg.norm(cellmat[0])
        b = np.linalg.norm(cellmat[1])
        c = np.linalg.norm(cellmat[2])
        ea = cellmat[0] / a
        eb = cellmat[1] / b
        ec = cellmat[2] / c
        A = np.degrees(np.arccos(eb @ ec))
        B = np.degrees(np.arccos(ec @ ea))
        C = np.degrees(np.arccos(ea @ eb))
        rotmat = np.linalg.inv(cellmat) @ cellvectors(a, b, c, A, B, C)
        logger.info("  The reshape matrix is reoriented.")
    else:
        rotmat = np.eye(3)

    atoms = []
    residue_count = 1
    for water in waters.values():
        for name, position in zip(water.labels, water.sites):
            atoms.append([residue_count, water.name, name, position])
        residue_count += 1
    for guest in guests:
        for name, position in zip(guest.labels, guest.sites):
            atoms.append([residue_count, guest.name, name, position])
        residue_count += 1
    for ion in ions.values():
        for name, position in zip(ion.labels, ion.sites):
            atoms.append([residue_count, ion.name, name, position])
        residue_count += 1

    logger.info("  Total number of atoms: {0}".format(len(atoms)))
    logger.info("  Total number of residues: {0}".format(residue_count))
    if len(atoms) > 99999:
        logger.warn(
            "  Fixed-digit format of Gromacs cannot deal with atoms more than 99999. Residue number and atom number are set appropriately."
        )
    s = ""
    s += "Generated by GenIce https://github.com/vitroid/GenIce \n"
    s += "{0}\n".format(len(atoms))
    for i, atom in enumerate(atoms):
        resno, resname, atomname, position = atom
        position = position @ rotmat
        if len(atoms) > 99999:
            s += "{0:5d}{1:5s}{2:>5s}{3:5d}{4:8.3f}{5:8.3f}{6:8.3f}\n".format(
                9999, resname, atomname, 9999, position[0], position[1], position[2]
            )
        else:
            s += "{0:5d}{1:5s}{2:>5s}{3:5d}{4:8.3f}{5:8.3f}{6:8.3f}\n".format(
                resno,
                resname,
                atomname,
                i + 1,
                position[0],
                position[1],
                position[2],
            )
    cellmat = cellmat @ rotmat
    if cellmat[1, 0] == 0 and cellmat[2, 0] == 0 and cellmat[2, 1] == 0:
        s += "    {0:.8f} {1:.8f} {2:.8f}\n".format(
            cellmat[0, 0], cellmat[1, 1], cellmat[2, 2]
        )
    else:
        s += "    {0:.8f} {1:.8f} {2:.8f} {3:.8f} {4:.8f} {5:.8f} {6:.8f} {7:.8f} {8:.8f}\n".format(
            cellmat[0, 0],
            cellmat[1, 1],
            cellmat[2, 2],
            cellmat[0, 1],
            cellmat[0, 2],
            cellmat[1, 0],
            cellmat[1, 2],
            cellmat[2, 0],
            cellmat[2, 1],
        )
    return s


if __name__ == "__main__":
    main()
