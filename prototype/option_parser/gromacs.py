"""
gromacs exporterプラグインのサンプル実装

このプラグインは guest, spot_guest, water_model オプションを処理し、
それ以外のオプションは返します。

注意: typeオプションはfoursite moleculeプラグインのオプションであり、
プラグインチェーンを通じてfoursiteプラグインに渡されます。
"""

from typing import Dict, Any, Tuple
from pool_parser import (
    parse_options_generic,
    OPTION_TYPE_STRING,
    OPTION_TYPE_KEYVALUE,
)


def parse_options(options: Dict[str, Any]) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    """
    gromacsプラグインのオプションを処理

    Args:
        options: オプションの辞書（設定ファイルの値が初期値として含まれる）

    Returns:
        (処理したオプション, 処理しなかったオプション)
    """
    # オプションの型定義
    # 注意: typeはfoursite moleculeプラグインのオプションなので、ここでは処理しない
    option_specs = {
        "guest": OPTION_TYPE_KEYVALUE,
        "spot_guest": OPTION_TYPE_KEYVALUE,
        "water_model": OPTION_TYPE_STRING,
    }

    return parse_options_generic(options, option_specs)


# テスト用
if __name__ == "__main__":
    # テストケース1: コマンドライン形式
    test_options1 = {
        "guest": ["A12=me", "A14=et"],
        "spot_guest": "0=foursite",
        "water_model": "foursite",
        "type": "ice",
        "shift": ("0.1", "0.1", "0.1"),  # これは処理しない
    }

    processed1, unprocessed1 = parse_options(test_options1)
    print("テストケース1:")
    print("処理したオプション:", processed1)
    print("処理しなかったオプション:", unprocessed1)
    print()

    # テストケース2: YAML形式（ネストされた構造）
    test_options2 = {
        "guest": {"A12": "me", "A14": "et"},
        "spot_guest": {"0": "foursite"},
        "water_model": {"name": "foursite", "options": {"type": "ice"}},
    }

    processed2, unprocessed2 = parse_options(test_options2)
    print("テストケース2:")
    print("処理したオプション:", processed2)
    print("処理しなかったオプション:", unprocessed2)
