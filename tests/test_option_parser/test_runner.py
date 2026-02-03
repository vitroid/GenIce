"""
プールベースパーサーのテストランナー（pytestなしで実行可能）
"""

import sys
import os

# パスを追加
sys.path.insert(0, os.path.dirname(__file__))

from pool_parser import PoolBasedParser
import tempfile


def test_parse_simple_args():
    """基本的な引数のパース"""
    print("Test: parse_simple_args")
    parser = PoolBasedParser()
    args = ["A15", "--rep", "2", "2", "2", "--seed", "42"]
    parser.parse_args(args)

    result = parser.get_result()
    assert (
        result["unitcell"]["name"] == "A15"
    ), f"Expected A15, got {result['unitcell']['name']}"
    assert result["base_options"]["replication_factors"] == ("2", "2", "2")
    assert result["base_options"]["seed"] == "42"
    print("  ✓ PASSED")


def test_parse_unitcell_options():
    """unitcellプラグインのオプションをパース"""
    print("Test: parse_unitcell_options")
    parser = PoolBasedParser()
    args = [
        "A15",
        "--rep",
        "2",
        "2",
        "2",
        "--shift",
        "0.1",
        "0.1",
        "0.1",
        "--density",
        "0.8",
        "--anion",
        "15=Cl",
        "--cation",
        "21=Na",
    ]
    parser.parse_args(args)

    result = parser.get_result()
    assert result["unitcell"]["name"] == "A15"
    assert result["unitcell"]["options"]["shift"] == ("0.1", "0.1", "0.1")
    assert result["unitcell"]["options"]["density"] == "0.8"
    assert result["unitcell"]["options"]["anion"] == "15=Cl"
    assert result["unitcell"]["options"]["cation"] == "21=Na"
    print("  ✓ PASSED")


def test_parse_complex_example():
    """複雑な例（OPTION_HANDLING_PLANS.mdの例）"""
    print("Test: parse_complex_example")
    parser = PoolBasedParser()
    args = [
        "A15",
        "--rep",
        "2",
        "2",
        "2",
        "--exporter",
        "gromacs",
        "--seed",
        "42",
        "--depol_loop",
        "2000",
        "--spot_anion",
        "1=Cl",
        "--spot_cation",
        "5=Na",
        "--shift",
        "0.1",
        "0.1",
        "0.1",
        "--anion",
        "15=Cl",
        "--cation",
        "21=Na",
        "--density",
        "0.8",
        "--guest",
        "A12=me",
        "--guest",
        "A14=et",
        "--spot_guest",
        "0=4site",
        "--water_model",
        "4site",
        "--type",
        "ice",
    ]
    parser.parse_args(args)

    result = parser.get_result()
    print(f"  unitcell name: {result['unitcell']['name']}")
    print(f"  base_options: {result['base_options']}")
    print(f"  unitcell_options: {result['unitcell']['options']}")
    print(f"  exporter name: {result['exporter']['name']}")

    assert result["unitcell"]["name"] == "A15"
    assert result["base_options"]["replication_factors"] == ("2", "2", "2")
    assert result["base_options"]["seed"] == "42"
    assert result["base_options"]["depol_loop"] == "2000"
    assert result["unitcell"]["options"]["shift"] == ("0.1", "0.1", "0.1")
    assert result["unitcell"]["options"]["density"] == "0.8"
    assert result["exporter"]["name"] == "gromacs"
    print("  ✓ PASSED")


def test_parse_yaml():
    """YAMLファイルからのパース"""
    print("Test: parse_yaml")
    try:
        import yaml
    except ImportError:
        print("  ⚠ SKIPPED (yaml module not available)")
        return

    yaml_content = """
genice3:
  seed: 42
  depol_loop: 2000
  replication_factors: [2, 2, 2]
  spot_anion:
    "1": Cl
  spot_cation:
    "5": Na

unitcell:
  name: A15
  shift: [0.1, 0.1, 0.1]
  anion:
    "15": Cl
  cation:
    "21": Na
  density: 0.8

exporter:
  name: gromacs
  guest:
    A12: me
    A14: et
  spot_guest:
    "0": "4site"
  water_model: "4site"
"""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
        f.write(yaml_content)
        yaml_path = f.name

    try:
        parser = PoolBasedParser()
        parser.parse_yaml(yaml_path)

        result = parser.get_result()
        assert result["base_options"]["seed"] == 42
        assert result["base_options"]["depol_loop"] == 2000
        assert result["base_options"]["replication_factors"] == (2, 2, 2)
        assert result["unitcell"]["name"] == "A15"
        assert result["unitcell"]["options"]["shift"] == [0.1, 0.1, 0.1]
        assert result["unitcell"]["options"]["density"] == 0.8
        assert result["exporter"]["name"] == "gromacs"
        assert result["exporter"]["options"]["guest"]["A12"] == "me"
        assert result["exporter"]["options"]["guest"]["A14"] == "et"
        print("  ✓ PASSED")
    finally:
        os.unlink(yaml_path)


def test_validate():
    """バリデーションテスト"""
    print("Test: validate")
    parser = PoolBasedParser()
    args = ["A15", "--rep", "2", "2", "2"]
    parser.parse_args(args)

    is_valid, errors = parser.validate()
    assert is_valid
    assert len(errors) == 0
    print("  ✓ PASSED")


def test_validate_missing_unitcell():
    """unitcell名がない場合のバリデーション"""
    print("Test: validate_missing_unitcell")
    parser = PoolBasedParser()
    args = ["--rep", "2", "2", "2"]
    parser.parse_args(args)

    is_valid, errors = parser.validate()
    assert not is_valid
    assert any("unitcell名" in error for error in errors)
    print("  ✓ PASSED")


if __name__ == "__main__":
    tests = [
        test_parse_simple_args,
        test_parse_unitcell_options,
        test_parse_complex_example,
        test_parse_yaml,
        test_validate,
        test_validate_missing_unitcell,
    ]

    passed = 0
    failed = 0

    for test in tests:
        try:
            test()
            passed += 1
        except Exception as e:
            print(f"  ✗ FAILED: {e}")
            import traceback

            traceback.print_exc()
            failed += 1

    print(f"\nResults: {passed} passed, {failed} failed")
