"""
A15 unitcellプラグインのサンプル実装

このプラグインは shift, density, anion, cation オプションを処理し、
それ以外のオプションは返します。
"""

from typing import Dict, Any, Tuple
from pool_parser import (
    parse_options_generic,
    OPTION_TYPE_STRING,
    OPTION_TYPE_TUPLE,
    OPTION_TYPE_KEYVALUE,
)


def parse_options(options: Dict[str, Any]) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    """
    A15プラグインのオプションを処理

    Args:
        options: オプションの辞書（設定ファイルの値が初期値として含まれる）

    Returns:
        (処理したオプション, 処理しなかったオプション)
    """
    # オプションの型定義
    option_specs = {
        "shift": OPTION_TYPE_TUPLE,
        "density": OPTION_TYPE_STRING,
        "anion": OPTION_TYPE_KEYVALUE,
        "cation": OPTION_TYPE_KEYVALUE,
    }

    # 後処理関数（数値への変換など）
    post_processors = {
        "shift": lambda x: [float(v) for v in x],
        "density": lambda x: float(x),
    }

    return parse_options_generic(options, option_specs, post_processors)


# テスト用
if __name__ == "__main__":
    # テストケース1: 基本的なケース
    test_options1 = {
        "shift": ("0.1", "0.1", "0.1"),
        "density": "0.8",
        "anion": "15=Cl",
        "cation": "21=Na",
        "guest": "A12=me",
        "water_model": "foursite",
    }

    processed1, unprocessed1 = parse_options(test_options1)
    print("テストケース1:")
    print("処理したオプション:", processed1)
    print("処理しなかったオプション:", unprocessed1)
    print()

    # テストケース2: 複数回指定されたanionとcation
    test_options2 = {
        "shift": ("0.1", "0.1", "0.1"),
        "density": "0.8",
        "anion": ["15=Cl", "16=F"],  # 複数回指定
        "cation": ["21=Na", "22=K"],  # 複数回指定
        "guest": "A12=me",
    }

    processed2, unprocessed2 = parse_options(test_options2)
    print("テストケース2（複数回指定）:")
    print("処理したオプション:", processed2)
    print("処理しなかったオプション:", unprocessed2)
