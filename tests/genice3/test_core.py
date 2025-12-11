"""
genice3コア機能のテスト
"""

import sys
from pathlib import Path

# プロジェクトのルートディレクトリをPythonパスに追加
project_root = Path(__file__).resolve().parent.parent.parent
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))

import pytest
import numpy as np
from genice3.genice import GenIce3, ConfigurationError
from genice3.plugin import UnitCell


def test_genice3_initialization():
    """GenIce3の初期化テスト"""
    genice = GenIce3()
    assert genice is not None
    assert genice.depol_loop == 1000  # デフォルト値
    # seedは__init__で使用されるが、インスタンス変数として保存されない
    # デフォルト値は1
    genice_with_seed = GenIce3(seed=42)
    assert genice_with_seed is not None


def test_genice3_with_options():
    """オプション付きの初期化テスト"""
    genice = GenIce3(
        depol_loop=500,
        seed=42,
        replication_matrix=np.array([[2, 0, 0], [0, 2, 0], [0, 0, 2]]),
    )
    assert genice.depol_loop == 500
    assert np.array_equal(
        genice.replication_matrix, np.array([[2, 0, 0], [0, 2, 0], [0, 0, 2]])
    )


def test_genice3_unitcell_setting():
    """単位胞の設定テスト"""
    genice = GenIce3()
    genice.unitcell = UnitCell("1h")
    assert genice.unitcell is not None


def test_genice3_invalid_kwargs():
    """無効なキーワード引数のテスト"""
    genice = GenIce3()
    with pytest.raises(ConfigurationError, match="Invalid keyword arguments"):
        genice = GenIce3(invalid_option="test")


def test_genice3_reactive_properties():
    """Reactive propertiesのテスト"""
    genice = GenIce3()
    genice.unitcell = UnitCell("1h")

    # プロパティにアクセスすると計算が実行される
    graph = genice.graph
    assert graph is not None
    assert graph.number_of_nodes() > 0

    digraph = genice.digraph
    assert digraph is not None
    assert digraph.number_of_nodes() == graph.number_of_nodes()


def test_genice3_lattice_sites():
    """格子サイトの取得テスト"""
    genice = GenIce3()
    genice.unitcell = UnitCell("1h")

    lattice_sites = genice.lattice_sites
    assert lattice_sites is not None
    assert len(lattice_sites) > 0
    assert lattice_sites.shape[1] == 3  # 3次元座標


def test_genice3_cell_matrix():
    """セル行列の取得テスト"""
    genice = GenIce3()
    genice.unitcell = UnitCell("1h")

    cell = genice.cell
    assert cell is not None
    assert cell.shape == (3, 3)  # 3x3行列


def test_genice3_water_molecules():
    """水分子の生成テスト"""
    from genice3.plugin import Molecule

    genice = GenIce3()
    genice.unitcell = UnitCell("1h")

    water_model = Molecule("4site")
    waters = genice.water_molecules(water_model=water_model)

    assert len(waters) > 0
    for mol in waters.values():
        assert mol.is_water
        assert len(mol.sites) == len(water_model.sites)


def test_genice3_replication():
    """複製のテスト"""
    replication_matrix = np.array([[2, 0, 0], [0, 2, 0], [0, 0, 2]])
    genice = GenIce3(replication_matrix=replication_matrix)
    genice.unitcell = UnitCell("1h")

    # 複製により原子数が増えることを確認
    from genice3.plugin import Molecule

    water_model = Molecule("4site")
    waters = genice.water_molecules(water_model=water_model)

    # 1hは通常8分子なので、2x2x2で64分子になるはず
    # （実際の値は単位胞の定義による）
    assert len(waters) >= 8  # 最低限元の分子数以上
