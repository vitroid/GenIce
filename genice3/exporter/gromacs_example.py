"""
Gromacs exporter pluginの例（デコレータ使用版）

このファイルは、`parse_plugin_options`デコレータの使用例を示すためのものです。
実際の実装では、このデコレータを使用することで、オプションの文字列パースが自動化されます。
"""

from logging import getLogger
from io import TextIOWrapper
import numpy as np
from typing import Dict, List, Any
from genice3.molecule import Molecule
from genice3.genice import GenIce3
from genice3.plugin import parse_plugin_options


def to_gro(
    cellmat: np.ndarray,
    waters: Dict[int, Molecule],
    guests: List[Molecule],
    ions: Dict[int, Molecule],
) -> str:
    """Gromacs形式で出力"""
    logger = getLogger("to_gro")
    logger.info("Generating .gro...")

    # ... (実装は省略) ...

    return ""


@parse_plugin_options
def dump(genice: GenIce3, file: TextIOWrapper, options: str = ""):
    """
    Gromacs形式で出力

    Args:
        genice: GenIce3インスタンス
        file: 出力先ファイル
        options: プラグインオプション（文字列形式）
            - 文字列で渡された場合、自動的に辞書にパースされる
            - 例: "guest.A12=me,shift=(0.1,0.1,0.1),verbose"
            - パース後: {"guest": {"A12": "me"}, "shift": [0.1, 0.1, 0.1], "verbose": True}

    使用例:
        # 呼び出し側
        dump(genice, sys.stdout, options="guest.A12=me,shift=(0.1,0.1,0.1)")

        # プラグイン側（自動的にパース済み）
        # options = {"guest": {"A12": "me"}, "shift": [0.1, 0.1, 0.1]}
    """
    # optionsは自動的に辞書に変換されている！
    # 文字列で渡された場合でも、ここでは既に辞書になっている

    # 例: ゲスト情報の取得
    guest_info = options.get("guest", {})
    if "A12" in guest_info:
        guest_spec = guest_info["A12"]
        # guest_specは文字列のまま（例: "0.3*me[monatomic]+0.6*et[molecular]"）
        # プラグイン側でさらに解釈が必要な場合は、ここで処理

    # 例: シフトベクトルの取得
    shift = options.get("shift", [0.0, 0.0, 0.0])
    # shiftは自動的にリストに変換されている（例: [0.1, 0.1, 0.1]）

    # 例: フラグの確認
    verbose = options.get("verbose", False)
    # verboseは自動的にTrueに変換されている

    # 構造データの取得
    structure = genice.get_atomic_structure()

    # 出力
    output = to_gro(
        cellmat=structure.cell,
        waters=structure.waters,
        guests=structure.guests,
        ions=structure.ions,
    )
    file.write(output)
