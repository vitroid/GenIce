"""
プールベースパーサーのテストランナー（pytestなしで実行可能）
"""

import sys
import os
import json

# パスを追加
sys.path.insert(0, os.path.dirname(__file__))

from pool_parser import PoolBasedParser
import tempfile

# プラグインをインポート（利用可能な場合）
try:
    from A15 import parse_options as parse_a15_options
    from gromacs import parse_options as parse_gromacs_options

    PLUGINS_AVAILABLE = True
except ImportError:
    PLUGINS_AVAILABLE = False

try:
    from foursite import parse_options as parse_foursite_options

    FOURSITE_AVAILABLE = True
except ImportError:
    FOURSITE_AVAILABLE = False

# グローバル変数: debugモード
DEBUG = False


def print_debug(message: str, indent: int = 0):
    """debugモードの時のみメッセージを表示"""
    if DEBUG:
        print("  " * indent + message)


def format_value(value):
    """値を読みやすい形式でフォーマット"""
    if isinstance(value, dict):
        return json.dumps(value, ensure_ascii=False, indent=2)
    elif isinstance(value, (list, tuple)):
        return json.dumps(list(value), ensure_ascii=False)
    else:
        return repr(value)


def test_parse_simple_args():
    """基本的な引数のパース"""
    print("Test: parse_simple_args")
    parser = PoolBasedParser()
    args = ["A15", "--rep", "2", "2", "2", "--seed", "42"]

    print_debug(f"入力引数: {args}")
    parser.parse_args(args)
    result = parser.get_result()

    print_debug("パース結果:", 1)
    print_debug(f"  unitcell名: {result['unitcell']['name']}", 2)
    print_debug(f"  基底レベルオプション:", 2)
    for key, value in result["base_options"].items():
        print_debug(f"    {key}: {format_value(value)}", 3)

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

    print_debug(f"入力引数: {args}")
    parser.parse_args(args)
    result = parser.get_result()

    print_debug("パース結果:", 1)
    print_debug(f"  unitcell名: {result['unitcell']['name']}", 2)
    print_debug(f"  unitcellオプション:", 2)
    for key, value in result["unitcell"]["options"].items():
        print_debug(f"    {key}: {format_value(value)} (型: {type(value).__name__})", 3)

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
        "0=foursite",
        "--water_model",
        "foursite",
        "--type",
        "ice",
    ]

    print_debug(f"入力引数: {args}")
    parser.parse_args(args)
    result = parser.get_result()

    print_debug("パース結果:", 1)
    print_debug(f"  unitcell名: {result['unitcell']['name']}", 2)
    print_debug(f"  exporter名: {result['exporter']['name']}", 2)
    print_debug(f"  基底レベルオプション:", 2)
    for key, value in result["base_options"].items():
        print_debug(f"    {key}: {format_value(value)}", 3)
    print_debug(f"  unitcellオプション:", 2)
    for key, value in result["unitcell"]["options"].items():
        print_debug(f"    {key}: {format_value(value)} (型: {type(value).__name__})", 3)
    print_debug(f"  exporterオプション:", 2)
    for key, value in result["exporter"]["options"].items():
        print_debug(f"    {key}: {format_value(value)}", 3)

    assert result["unitcell"]["name"] == "A15"
    assert result["base_options"]["replication_factors"] == ("2", "2", "2")
    assert result["base_options"]["seed"] == "42"
    assert result["base_options"]["depol_loop"] == "2000"
    assert result["unitcell"]["options"]["shift"] == ("0.1", "0.1", "0.1")
    assert result["unitcell"]["options"]["density"] == "0.8"
    assert result["exporter"]["name"] == "gromacs"

    # gromacsプラグインのオプションがunitcell_optionsに含まれていることを確認
    # （unitcellプラグインが処理しなかったものはexporterに渡される）
    assert (
        "guest" in result["unitcell"]["options"]
    ), "guestオプションがunitcell_optionsに含まれる"
    assert (
        "spot_guest" in result["unitcell"]["options"]
    ), "spot_guestオプションがunitcell_optionsに含まれる"
    assert (
        "water_model" in result["unitcell"]["options"]
    ), "water_modelオプションがunitcell_optionsに含まれる"
    assert (
        "type" in result["unitcell"]["options"]
    ), "typeオプションがunitcell_optionsに含まれる"

    # プラグインが利用可能な場合、実際にプラグインで処理して確認
    if PLUGINS_AVAILABLE:
        print("\n  A15プラグインでの処理:")
        print(f"    入力オプション: {list(result['unitcell']['options'].keys())}")

        # A15プラグインで処理
        a15_processed, a15_unprocessed = parse_a15_options(
            result["unitcell"]["options"]
        )

        print(f"    処理したオプション:")
        for key, value in a15_processed.items():
            original_value = result["unitcell"]["options"].get(key, "N/A")
            print(f"      {key}: {original_value} → {format_value(value)}")

        print(f"    処理しなかったオプション（exporterに渡す）:")
        for key, value in a15_unprocessed.items():
            print(f"      {key}: {value}")

        # gromacsプラグインのオプションがA15で処理されず、exporterに渡されることを確認
        assert (
            "guest" in a15_unprocessed
        ), "guestオプションはA15で処理されず、gromacsに渡される"
        assert (
            "spot_guest" in a15_unprocessed
        ), "spot_guestオプションはA15で処理されず、gromacsに渡される"
        assert (
            "water_model" in a15_unprocessed
        ), "water_modelオプションはA15で処理されず、gromacsに渡される"
        assert (
            "type" in a15_unprocessed
        ), "typeオプションはA15で処理されず、gromacsに渡される"

        # gromacsプラグインで処理
        gromacs_options = {**result["exporter"]["options"], **a15_unprocessed}
        print(f"\n  gromacsプラグインでの処理:")
        print(f"    入力オプション: {list(gromacs_options.keys())}")

        gromacs_processed, gromacs_unprocessed = parse_gromacs_options(gromacs_options)

        print(f"    処理したオプション:")
        for key, value in gromacs_processed.items():
            original_value = gromacs_options.get(key, "N/A")
            print(
                f"      {key}: {format_value(original_value)} → {format_value(value)}"
            )

        if gromacs_unprocessed:
            print(f"    処理しなかったオプション:")
            for key, value in gromacs_unprocessed.items():
                print(f"      {key}: {value}")
        else:
            print(f"    処理しなかったオプション: なし")

        # foursiteプラグインをチェーンとして処理
        # water_modelがfoursiteの場合のみ、foursiteプラグインを呼び出す
        foursite_processed = {}
        if (
            FOURSITE_AVAILABLE
            and gromacs_unprocessed
            and "type" in gromacs_unprocessed
            and "water_model" in gromacs_processed
        ):
            water_model_name = gromacs_processed["water_model"]
            if isinstance(water_model_name, dict):
                water_model_name = water_model_name.get("name")
            if water_model_name == "foursite":
                # gromacsで処理したオプション（water_model）も含めてfoursiteに渡す
                foursite_options = {**gromacs_processed, **gromacs_unprocessed}
                print(f"\n  foursite moleculeプラグインでの処理（チェーン）:")
                print(f"    入力オプション: {list(foursite_options.keys())}")
                foursite_processed, foursite_unprocessed = parse_foursite_options(
                    foursite_options
                )

                print(f"    処理したオプション:")
                for key, value in foursite_processed.items():
                    original_value = foursite_options.get(key, "N/A")
                    print(
                        f"      {key}: {format_value(original_value)} → {format_value(value)}"
                    )

                if foursite_unprocessed:
                    print(f"    処理しなかったオプション:")
                    for key, value in foursite_unprocessed.items():
                        print(f"      {key}: {value}")
                else:
                    print(f"    処理しなかったオプション: なし")

                # foursiteで処理したオプションをgromacsの処理結果に統合
                if "water_model" in gromacs_processed and "type" in foursite_processed:
                    if isinstance(gromacs_processed["water_model"], str):
                        gromacs_processed["water_model"] = {
                            "name": gromacs_processed["water_model"],
                            "options": foursite_processed,
                        }
                    elif isinstance(gromacs_processed["water_model"], dict):
                        if "options" not in gromacs_processed["water_model"]:
                            gromacs_processed["water_model"]["options"] = {}
                        gromacs_processed["water_model"]["options"].update(
                            foursite_processed
                        )

        # gromacsプラグインがこれらのオプションを処理することを確認
        assert (
            "guest" in gromacs_processed
        ), "guestオプションはgromacsプラグインで処理される"
        assert (
            "spot_guest" in gromacs_processed
        ), "spot_guestオプションはgromacsプラグインで処理される"
        assert (
            "water_model" in gromacs_processed
        ), "water_modelオプションはgromacsプラグインで処理される"
        # typeはfoursite moleculeプラグインのオプションで、water_modelがfoursiteの場合に
        # foursiteプラグインが処理するため、water_modelの中に含まれる
        assert isinstance(
            gromacs_processed["water_model"], dict
        ), "water_modelは辞書として処理される"
        assert (
            "options" in gromacs_processed["water_model"]
        ), "water_modelにoptionsが含まれる"
        # foursiteプラグインがチェーンとして処理されたことを確認
        assert (
            "type" in gromacs_processed["water_model"]["options"]
        ), "typeはfoursite moleculeプラグインのオプションとして、water_modelのoptionsの中に含まれる"

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
    # コマンドライン引数からdebugモードを判定
    if "--debug" in sys.argv or "-d" in sys.argv:
        DEBUG = True
        print("Debug mode: ON")
        print()
    else:
        DEBUG = False

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
