"""
foursite moleculeプラグインのサンプル実装

このプラグインは type オプションを処理します。
自分が呼び出されたということは、water_modelがfoursiteであることが決定されているため、
water_modelの値を確認する必要はありません。
"""

from typing import Dict, Any, Tuple
from pool_parser import parse_options_generic, OPTION_TYPE_STRING


def parse_options(options: Dict[str, Any]) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    """
    foursiteプラグインのオプションを処理

    Args:
        options: オプションの辞書（typeを含む）

    Returns:
        (処理したオプション, 処理しなかったオプション)
    """
    # オプションの型定義
    option_specs = {
        "type": OPTION_TYPE_STRING,
    }

    return parse_options_generic(options, option_specs)


# テスト用
if __name__ == "__main__":
    # テストケース1: typeオプションを処理
    test_options1 = {
        "type": "ice",
        "other": "value",  # これは処理しない
    }

    processed1, unprocessed1 = parse_options(test_options1)
    print("テストケース1:")
    print("処理したオプション:", processed1)
    print("処理しなかったオプション:", unprocessed1)
    print()

    # テストケース2: typeオプションを処理（別の値）
    test_options2 = {
        "type": "liquid",
    }

    processed2, unprocessed2 = parse_options(test_options2)
    print("テストケース2:")
    print("処理したオプション:", processed2)
    print("処理しなかったオプション:", unprocessed2)
