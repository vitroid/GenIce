"""
Pytest fixtures for genice3 tests
"""

import sys
from pathlib import Path
import os

# プロジェクトのルートディレクトリをPythonパスに追加
# __file__の絶対パスから、tests/genice3/conftest.py -> tests/genice3 -> tests -> . と遡る
try:
    # __file__が定義されている場合（通常の実行）
    project_root = Path(__file__).resolve().parent.parent.parent
except NameError:
    # __file__が定義されていない場合（exec等で実行された場合）
    # カレントディレクトリから推測
    project_root = Path.cwd()

project_root_str = str(project_root)
if project_root_str not in sys.path:
    sys.path.insert(0, project_root_str)

# 環境変数PYTHONPATHにも設定（確実にするため）
if "PYTHONPATH" not in os.environ:
    os.environ["PYTHONPATH"] = project_root_str
elif project_root_str not in os.environ["PYTHONPATH"].split(os.pathsep):
    os.environ["PYTHONPATH"] = project_root_str + os.pathsep + os.environ["PYTHONPATH"]

import pytest
import numpy as np
from io import StringIO


# 遅延インポート: fixture内でインポートする
def _get_genice3():
    """GenIce3クラスを取得"""
    from genice3.genice import GenIce3

    return GenIce3


def _get_unitcell():
    """UnitCell関数を取得"""
    from genice3.plugin import UnitCell

    return UnitCell


def _get_exporter():
    """Exporter関数を取得"""
    from genice3.plugin import Exporter

    return Exporter


def _get_molecule():
    """Molecule関数を取得"""
    from genice3.plugin import Molecule

    return Molecule


@pytest.fixture
def basic_genice():
    """基本的なGenIce3インスタンス"""
    GenIce3 = _get_genice3()
    return GenIce3(seed=1)


@pytest.fixture
def genice_1h(basic_genice):
    """1h単位胞を持つGenIce3インスタンス"""
    UnitCell = _get_unitcell()
    basic_genice.unitcell = UnitCell("1h")
    return basic_genice


@pytest.fixture
def genice_A15(basic_genice):
    """A15単位胞を持つGenIce3インスタンス"""
    UnitCell = _get_unitcell()
    basic_genice.unitcell = UnitCell("A15")
    return basic_genice


@pytest.fixture
def water_4site():
    """4site水分子モデル"""
    Molecule = _get_molecule()
    return Molecule("4site")


@pytest.fixture
def water_tip4p():
    """TIP4P水分子モデル"""
    Molecule = _get_molecule()
    return Molecule("tip4p")


@pytest.fixture
def exporter_gromacs():
    """Gromacs exporter"""
    Exporter = _get_exporter()
    return Exporter("gromacs")


@pytest.fixture
def exporter_cif():
    """CIF exporter"""
    Exporter = _get_exporter()
    return Exporter("cif")


@pytest.fixture
def output_file():
    """出力用のStringIO"""
    return StringIO()
