"""
YAML設定ファイルの読み込みと処理

設定ファイルからオプションを読み込み、コマンドライン引数と統合します。
コマンドライン引数が設定ファイルの値を上書きします。
"""

import sys
from pathlib import Path
from typing import Dict, Any, Optional
from logging import getLogger
import numpy as np

logger = getLogger(__name__)

# YAMLライブラリのインポート
try:
    import yaml

    YAML_AVAILABLE = True
except ImportError:
    YAML_AVAILABLE = False


def load_config(config_path: str) -> Dict[str, Any]:
    """
    YAML設定ファイルを読み込む

    Args:
        config_path: 設定ファイルのパス

    Returns:
        設定の辞書。設定ファイルが存在しない、または読み込めない場合は空の辞書を返す。

    Raises:
        FileNotFoundError: 設定ファイルが存在しない場合
        ValueError: YAMLライブラリが利用できない場合
    """
    if not YAML_AVAILABLE:
        raise ValueError(
            "YAMLライブラリが利用できません。`pip install pyyaml`でインストールしてください。"
        )

    config_file = Path(config_path)
    if not config_file.exists():
        raise FileNotFoundError(f"設定ファイルが見つかりません: {config_path}")

    logger.debug(f"設定ファイルを読み込み中: {config_path}")

    try:
        with open(config_file, "r", encoding="utf-8") as f:
            config = yaml.safe_load(f)
        logger.debug(f"設定ファイルを正常に読み込みました: {config}")
        return config or {}
    except Exception as e:
        logger.error(f"設定ファイルの読み込みに失敗しました: {e}")
        raise


def parse_config(config: Dict[str, Any]) -> Dict[str, Any]:
    """
    YAML設定ファイルの内容をGenIce3のオプション形式に変換

    Args:
        config: YAML設定ファイルから読み込んだ辞書

    Returns:
        GenIce3のオプション形式の辞書

    設定ファイルの構造例:
        genice3:
          seed: 1
          depol_loop: 1000
          replication_factors: [2, 2, 2]
          assess_cages: false
            spot_anion:
              "1": "Cl"

            spot_cation:
              "1": "Na"

        unitcell:
          name: "1h"
          shift: [0.1, 0.1, 0.1]
          # その他のオプションはプラグイン側で解釈

        exporter:
          name: "gromacs"
          guest:
            A12: "me"
          # その他のオプションはプラグイン側で解釈

    """
    result: Dict[str, Any] = {}

    # genice3セクション（メインオプション）
    if "genice3" in config:
        genice3_config = config["genice3"]
        if "seed" in genice3_config:
            result["seed"] = genice3_config["seed"]
        if "depol_loop" in genice3_config:
            result["depol_loop"] = genice3_config["depol_loop"]
        if "replication_factors" in genice3_config:
            result["replication_factors"] = tuple(genice3_config["replication_factors"])
        if "replication_matrix" in genice3_config:
            matrix = genice3_config["replication_matrix"]
            if isinstance(matrix, list) and len(matrix) == 9:
                result["replication_matrix"] = np.array(matrix).reshape(3, 3)
        if "assess_cages" in genice3_config:
            result["assess_cages"] = genice3_config["assess_cages"]
        if "debug" in genice3_config:
            result["debug"] = genice3_config["debug"]
        if "spot_anion" in genice3_config:
            result["spot_anion"] = genice3_config["spot_anion"]
        if "spot_cation" in genice3_config:
            result["spot_cation"] = genice3_config["spot_cation"]

    # unitcellセクションとexporterセクションはそのまま返す（プラグイン側で処理）
    if "unitcell" in config:
        result["unitcell"] = config["unitcell"]
    if "exporter" in config:
        result["exporter"] = config["exporter"]

    return result


def _is_int_key(key: Any) -> bool:
    """
    キーが整数に変換可能かどうかを判定

    Args:
        key: チェックするキー

    Returns:
        整数に変換可能な場合True
    """
    try:
        int(str(key))
        return True
    except (ValueError, TypeError):
        return False


def merge_config_with_args(
    config_options: Dict[str, Any], args: Dict[str, Any]
) -> Dict[str, Any]:
    """
    設定ファイルのオプションとコマンドライン引数を統合
    コマンドライン引数が設定ファイルの値を上書きします。

    Args:
        config_options: 設定ファイルから読み込んだオプション
        args: コマンドライン引数

    Returns:
        統合されたオプション
    """
    merged = config_options.copy()

    logger.info(f"config: {config_options}")
    logger.info(f"option: {args}")

    # コマンドライン引数で上書き
    for key, value in args.items():
        if value is not None:
            # リストの場合は結合（spot_anion, spot_cationなど）
            if isinstance(value, dict) and len(value) > 0:
                merged[key] = value
            else:
                merged[key] = value
    logger.info(f"merged: {merged}")
    return merged
