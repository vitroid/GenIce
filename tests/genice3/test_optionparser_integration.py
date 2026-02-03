"""
genice3/pool_parser.pyの統合テスト

PoolBasedParserの統合が正しく動作するかを確認するテスト
"""

import sys
from pathlib import Path

# プロジェクトルートをパスに追加
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

from genice3.cli.pool_parser import PoolBasedParser


def test_basic_parsing():
    """基本的なパースのテスト"""
    parser = PoolBasedParser()
    parser.parse_args(["A15", "--exporter", "gromacs", "--rep", "2", "2", "2"])
    result = parser.get_result()

    assert result["unitcell"]["name"] == "A15"
    assert result["exporter"]["name"] == "gromacs"
    assert result["base_options"]["replication_factors"] == ("2", "2", "2")
    print("✓ 基本パーステスト成功")


def test_complex_parsing():
    """複雑なパースのテスト"""
    parser = PoolBasedParser()
    parser.parse_args(
        [
            "A15",
            "--exporter",
            "gromacs",
            "--rep",
            "2",
            "2",
            "2",
            "--seed",
            "42",
            "--spot_anion",
            "1=Cl",
            "--spot_cation",
            "5=Na",
        ]
    )
    result = parser.get_result()

    assert result["unitcell"]["name"] == "A15"
    assert result["exporter"]["name"] == "gromacs"
    assert result["base_options"]["seed"] == "42"
    assert result["base_options"]["spot_anion"] == {"1": "Cl"}
    assert result["base_options"]["spot_cation"] == {"5": "Na"}
    print("✓ 複雑なパーステスト成功")


def test_validation():
    """バリデーションのテスト"""
    parser = PoolBasedParser()
    parser.parse_args(["A15", "--exporter", "gromacs"])
    is_valid, errors = parser.validate()
    assert is_valid, f"バリデーションエラー: {errors}"
    print("✓ バリデーションテスト成功")


def test_missing_unitcell():
    """unitcellが指定されていない場合のテスト"""
    parser = PoolBasedParser()
    parser.parse_args(["--exporter", "gromacs"])
    is_valid, errors = parser.validate()
    assert not is_valid, "unitcellが指定されていない場合はエラーになるべき"
    assert any("unitcell" in error.lower() for error in errors)
    print("✓ unitcell未指定のテスト成功")


if __name__ == "__main__":
    print("=" * 60)
    print("genice3/pool_parser.py 統合テスト")
    print("=" * 60)
    print()

    try:
        test_basic_parsing()
        test_complex_parsing()
        test_validation()
        test_missing_unitcell()
        print()
        print("=" * 60)
        print("✓ すべてのテストが成功しました")
        print("=" * 60)
    except Exception as e:
        print(f"✗ テスト失敗: {e}")
        import traceback

        traceback.print_exc()
        sys.exit(1)
