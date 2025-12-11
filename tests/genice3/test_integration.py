"""
genice3統合テスト
代表的な組み合わせで動作確認
"""

import sys
from pathlib import Path

# プロジェクトのルートディレクトリをPythonパスに追加
project_root = Path(__file__).resolve().parent.parent.parent
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))

import pytest
from io import StringIO
from genice3.genice import GenIce3
from genice3.plugin import UnitCell, Exporter, Molecule
from tests.genice3.test_validation import (
    validate_comprehensive,
    validate_basic_execution,
    validate_output_format,
    validate_ice_rules,
    validate_atom_count,
    validate_coordinates_in_cell,
    validate_graph_structure,
    validate_reproducibility,
)


# テストマトリックス: 代表的な組み合わせ
TEST_COMBINATIONS = [
    ("1h", "4site", "gromacs"),
    ("1h", "tip4p", "gromacs"),
    ("A15", "4site", "gromacs"),
    ("1h", "4site", "cif"),
    ("A15", "4site", "cif"),
]


@pytest.mark.parametrize("unitcell_name,molecule_name,exporter_name", TEST_COMBINATIONS)
def test_basic_combination(unitcell_name, molecule_name, exporter_name):
    """基本的な組み合わせが動作することを確認（レベル1+2）"""
    genice = GenIce3(seed=1)
    genice.unitcell = UnitCell(unitcell_name)

    exporter_module = Exporter(exporter_name)
    water_model = Molecule(molecule_name)

    output_file = StringIO()
    # water_model_processorはプラグイン名（例: "4site"）を期待するため、
    # water_model.name（例: "ICE"）ではなく、元のプラグイン名を渡す
    exporter_module.dump(genice, output_file, water_model=molecule_name)
    output = output_file.getvalue()

    # レベル1: 基本動作
    validate_basic_execution(genice, exporter_module, output)

    # レベル2: 出力形式
    validate_output_format(exporter_name, output)


@pytest.mark.parametrize("unitcell_name", ["1h", "A15"])
@pytest.mark.parametrize("molecule_name", ["4site", "tip4p"])
def test_ice_rules_validation(unitcell_name, molecule_name):
    """アイスルールの検証（レベル3）"""
    genice = GenIce3(seed=1)
    genice.unitcell = UnitCell(unitcell_name)

    water_model = Molecule(molecule_name)

    # レベル3: アイスルール
    validate_ice_rules(genice)

    # 原子数と座標の検証
    validate_atom_count(genice, water_model)
    validate_coordinates_in_cell(genice, water_model)


@pytest.mark.parametrize("unitcell_name", ["1h", "A15"])
def test_graph_structure_validation(unitcell_name):
    """グラフ構造の検証（レベル5）"""
    genice = GenIce3(seed=1)
    genice.unitcell = UnitCell(unitcell_name)

    validate_graph_structure(genice)


@pytest.mark.parametrize("unitcell_name", ["1h", "A15"])
def test_reproducibility(unitcell_name):
    """再現性のテスト（レベル4）"""
    validate_reproducibility(unitcell_name, seed=1)


@pytest.mark.parametrize(
    "unitcell_name,molecule_name,exporter_name",
    [("1h", "4site", "gromacs"), ("A15", "4site", "cif")],
)
def test_comprehensive_validation(unitcell_name, molecule_name, exporter_name):
    """総合的な検証"""
    genice = GenIce3(seed=1)
    genice.unitcell = UnitCell(unitcell_name)

    exporter_module = Exporter(exporter_name)
    water_model = Molecule(molecule_name)

    output_file = StringIO()
    # water_model_processorはプラグイン名（例: "4site"）を期待するため、
    # water_model.name（例: "ICE"）ではなく、元のプラグイン名を渡す
    exporter_module.dump(genice, output_file, water_model=molecule_name)
    output = output_file.getvalue()

    # 総合検証
    validate_comprehensive(genice, exporter_name, output, water_model)


def test_exporter_options():
    """Exporterオプションのテスト"""
    genice = GenIce3(seed=1)
    genice.unitcell = UnitCell("A15")

    exporter_module = Exporter("gromacs")
    output_file = StringIO()

    # オプション付きで実行
    exporter_module.dump(
        genice,
        output_file,
        water_model="4site",
        guest={"A12": "me"},
        spot_guest={0: "4site"},
    )

    output = output_file.getvalue()
    assert len(output) > 0


def test_unitcell_options():
    """UnitCellオプションのテスト"""
    genice = GenIce3(seed=1)

    # オプション付きで単位胞を設定
    genice.unitcell = UnitCell("A15", shift=(0.1, 0.1, 0.1), density=0.8)

    assert genice.unitcell is not None

    # グラフが生成できることを確認
    graph = genice.graph
    assert graph.number_of_nodes() > 0


def test_ion_doping():
    """イオン置換のテスト"""
    genice = GenIce3(seed=1)
    genice.unitcell = UnitCell("1h")

    # イオンを置換
    genice.spot_anions = {0: "Cl"}
    genice.spot_cations = {1: "Na"}

    # グラフが生成できることを確認
    graph = genice.graph
    assert graph.number_of_nodes() > 0
