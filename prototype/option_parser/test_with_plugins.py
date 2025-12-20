"""
プールベースオプションパーサーのテスト（プラグインを使用）

実際のプラグインを使用して、オプションが正しく処理されることを確認します。
"""

import sys
import os
import json

# パスを追加
sys.path.insert(0, os.path.dirname(__file__))

from pool_parser import PoolBasedParser
from A15 import parse_options as parse_a15_options
from gromacs import parse_options as parse_gromacs_options

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


def test_parse_with_plugins():
    """プラグインを使用したパーステスト"""
    print("=" * 60)
    print("プラグインを使用したパーステスト")
    print("=" * 60)

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

    print("\n基底レベルのオプション:")
    for key, value in result["base_options"].items():
        print(f"  {key}: {value}")
        print_debug(f"    → 型: {type(value).__name__}, 値: {format_value(value)}", 2)

    print(f"\nunitcell: {result['unitcell']['name']}")
    print("unitcellに渡されるオプション:")
    print(f"  {result['unitcell']['options']}")
    print_debug("unitcellオプションの詳細:", 1)
    for key, value in result["unitcell"]["options"].items():
        print_debug(f"  {key}: {format_value(value)} (型: {type(value).__name__})", 2)

    # A15プラグインでオプションを処理
    print_debug("\nA15プラグインでの処理:", 0)
    print_debug(f"  入力オプション: {format_value(result['unitcell']['options'])}", 1)
    a15_processed, a15_unprocessed = parse_a15_options(result["unitcell"]["options"])

    print("\nA15プラグインが処理したオプション:")
    for key, value in a15_processed.items():
        print(f"  {key}: {value}")
        print_debug(f"    → 元の値: {result['unitcell']['options'].get(key, 'N/A')}", 3)
        print_debug(
            f"    → 変換後: {format_value(value)} (型: {type(value).__name__})", 3
        )

    print("\nA15プラグインが処理しなかったオプション（exporterに渡す）:")
    for key, value in a15_unprocessed.items():
        print(f"  {key}: {value}")

    print(f"\nexporter: {result['exporter']['name']}")
    print("exporterに渡されるオプション（A15で処理されなかったもの）:")
    print(f"  {a15_unprocessed}")

    # gromacsプラグインでオプションを処理
    # まず、元のexporterオプションと統合
    gromacs_options = {**result["exporter"]["options"], **a15_unprocessed}
    print_debug("\ngromacsプラグインでの処理:", 0)
    print_debug(f"  入力オプション: {format_value(gromacs_options)}", 1)
    gromacs_processed, gromacs_unprocessed = parse_gromacs_options(gromacs_options)

    print("\ngromacsプラグインが処理したオプション:")
    for key, value in gromacs_processed.items():
        print(f"  {key}: {value}")
        print_debug(f"    → 元の値: {gromacs_options.get(key, 'N/A')}", 3)
        print_debug(
            f"    → 変換後: {format_value(value)} (型: {type(value).__name__})", 3
        )

    print("\ngromacsプラグインが処理しなかったオプション:")
    for key, value in gromacs_unprocessed.items():
        print(f"  {key}: {value}")

    # foursiteプラグインをチェーンとして処理（gromacsが処理しなかったオプションを渡す）
    try:
        from foursite import parse_options as parse_foursite_options

        FOURSITE_AVAILABLE = True
    except ImportError:
        FOURSITE_AVAILABLE = False

    # foursiteプラグインをチェーンとして処理
    # water_modelがfoursiteの場合のみ、foursiteプラグインを呼び出す
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
            print_debug("\nfoursiteプラグインでの処理:", 0)
            print_debug(f"  入力オプション: {format_value(foursite_options)}", 1)
            foursite_processed, foursite_unprocessed = parse_foursite_options(
                foursite_options
            )

            if foursite_processed:
                print("\nfoursiteプラグインが処理したオプション:")
                for key, value in foursite_processed.items():
                    print(f"  {key}: {value}")
                    print_debug(
                        f"    → 変換後: {format_value(value)} (型: {type(value).__name__})",
                        3,
                    )
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

            if foursite_unprocessed:
            print("\nfoursiteプラグインが処理しなかったオプション:")
            for key, value in foursite_unprocessed.items():
                print(f"  {key}: {value}")
            # foursiteで処理しなかったものは、最終的な未処理オプション
            final_unprocessed = foursite_unprocessed
        else:
            final_unprocessed = {}
    else:
        final_unprocessed = gromacs_unprocessed

    # 最終的に処理されなかったオプションがあればエラー
    if final_unprocessed:
        print("\n⚠ 警告: 処理されなかったオプションがあります:")
        for key, value in final_unprocessed.items():
            print(f"  {key}: {value}")
    else:
        print("\n✓ すべてのオプションが処理されました")

    # アサーション
    assert result["unitcell"]["name"] == "A15"
    assert "shift" in a15_processed
    assert "density" in a15_processed
    assert "anion" in a15_processed
    assert "cation" in a15_processed

    assert "guest" in gromacs_processed or "guest" in a15_unprocessed
    assert "water_model" in gromacs_processed or "water_model" in a15_unprocessed

    print("\n✓ テスト成功")


def test_parse_simple_with_plugins():
    """シンプルなケースのテスト"""
    print("\n" + "=" * 60)
    print("シンプルなケースのテスト")
    print("=" * 60)

    parser = PoolBasedParser()
    args = [
        "A15",
        "--shift",
        "0.1",
        "0.1",
        "0.1",
        "--density",
        "0.8",
        "--guest",
        "A12=me",
    ]
    parser.parse_args(args)

    result = parser.get_result()

    # A15プラグインで処理
    a15_processed, a15_unprocessed = parse_a15_options(result["unitcell"]["options"])

    print("\nA15プラグインが処理したオプション:")
    for key, value in a15_processed.items():
        print(f"  {key}: {value}")

    print("\nA15プラグインが処理しなかったオプション:")
    for key, value in a15_unprocessed.items():
        print(f"  {key}: {value}")

    # guestオプションはA15では処理されないはず
    assert "guest" in a15_unprocessed

    # gromacsプラグインで処理（exporterが指定されていない場合はスキップ）
    if result["exporter"]["name"]:
        gromacs_options = {**result["exporter"]["options"], **a15_unprocessed}
        gromacs_processed, gromacs_unprocessed = parse_gromacs_options(gromacs_options)
        print("\ngromacsプラグインが処理したオプション:")
        for key, value in gromacs_processed.items():
            print(f"  {key}: {value}")

    print("\n✓ テスト成功")


if __name__ == "__main__":
    # コマンドライン引数からdebugモードを判定
    if "--debug" in sys.argv or "-d" in sys.argv:
        DEBUG = True
        print("Debug mode: ON")
        print()
    else:
        DEBUG = False

    test_parse_with_plugins()
    test_parse_simple_with_plugins()
