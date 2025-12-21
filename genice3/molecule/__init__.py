from dataclasses import dataclass
import numpy as np
from math import sin, cos, radians
from typing import Dict, Any, Tuple
from genice3.cli.pool_parser import (
    parse_options_generic,
    OPTION_TYPE_STRING,
)


@dataclass
class Molecule:
    """
    Base class of a molecule
    """

    sites: np.ndarray
    labels: list[str]
    name: str
    is_water: bool = False

    def __repr__(self) -> str:
        return (
            f"Molecule(name={self.name!r}, "
            f"n_sites={len(self.sites)}, "
            f"labels={self.labels}, "
            f"is_water={self.is_water})"
        )

    @staticmethod
    def parse_options(options: Dict[str, Any]) -> Tuple[Dict[str, Any], Dict[str, Any]]:
        """
        moleculeプラグインの共通オプションを処理

        この関数は、動的プラグインチェーン実行システムから呼び出されます。
        すべてのmoleculeプラグインで共通のオプションを処理します。

        現在、共通オプションは特にありませんが、将来的に共通オプションが追加された場合に備えて
        このメソッドを提供しています。個別のmoleculeプラグインが独自のオプションを持つ場合は、
        プラグイン固有のparse_options関数を定義することで、基底クラスのメソッドを上書きできます。

        Args:
            options: オプションの辞書

        Returns:
            (処理したオプション, 処理しなかったオプション) のタプル
            - 処理したオプション: moleculeプラグインが処理したオプション（現在は空）
            - 処理しなかったオプション: その他のオプション
        """
        # 現在、moleculeプラグインで共通のオプションはないため、空の辞書を返す
        # 個別のプラグインが独自のオプションを持つ場合は、プラグイン固有のparse_options関数を定義
        return ({}, options)
