from logging import getLogger, basicConfig, DEBUG, INFO
from typing import Optional, Tuple, Dict, Any
import numpy as np
import sys
from typing import Dict, List
from importlib.metadata import version, PackageNotFoundError

from genice3.plugin import safe_import
from genice3.genice import GenIce3
from genice3.cli.pool_parser import (
    PoolBasedParser,
    parse_options_generic,
    OPTION_TYPE_FLAG,
    OPTION_TYPE_STRING,
    OPTION_TYPE_TUPLE,
    OPTION_TYPE_KEYVALUE,
)
from genice3.unitcell import ion_processor

HELP_VERSION = "Show the version and exit."


def get_version():
    """パッケージのバージョンを取得する。見つからない場合はデフォルト値を返す。"""
    try:
        return version("genice3")
    except PackageNotFoundError:
        return "3.*.*.*"


def parse_base_options(
    options: Dict[str, Any],
) -> Dict[str, Any]:
    """
    GenIce3のbaseレベルオプションを処理

    この関数は、動的プラグインチェーン実行システムから呼び出されます。
    GenIce3クラスが受け取るbaseレベルのオプション（seed, depol_loop, replication_matrix,
    replication_factors, assess_cages, debug, spot_anion, spot_cation）を処理します。

    Args:
        options: オプションの辞書（設定ファイルの値が初期値として含まれる可能性がある）
            - seed: 乱数シード（文字列形式、例: "123"）
            - depol_loop: 双極子最適化の反復回数（文字列形式、例: "1000"）
            - replication_factors: 複製ファクター（タプル形式、例: (2, 2, 2)）
            - replication_matrix: 複製行列（タプル形式、例: (2, 0, 0, 0, 2, 0, 0, 0, 2)）
            - assess_cages: ケージ調査を行うかどうか（文字列形式、例: "true"）
            - debug: デバッグモード（文字列形式、例: "true"）
            - spot_anion: アニオン配置（key=value形式、例: "15=Cl" または {"15": "Cl"}）
            - spot_cation: カチオン配置（key=value形式、例: "21=Na" または {"21": "Na"}）

    Returns:
        統合されたオプション辞書（baseレベルで処理されたオプション（型変換済み）と
        処理されなかったオプション（プラグインに渡す）の両方を含む）
    """
    option_specs = {
        "seed": OPTION_TYPE_STRING,  # 文字列として取得し、後処理で整数に変換
        "depol_loop": OPTION_TYPE_STRING,  # 文字列として取得し、後処理で整数に変換
        "replication_factors": OPTION_TYPE_TUPLE,  # タプルとして取得し、後処理で整数タプルに変換
        "replication_matrix": OPTION_TYPE_TUPLE,  # タプルとして処理（9要素のリスト）
        "assess_cages": OPTION_TYPE_STRING,  # 文字列として取得し、後処理でブール値に変換（Falseの値を正しく処理するため）
        "debug": OPTION_TYPE_FLAG,  # フラグ型（引数なし）
        "spot_anion": OPTION_TYPE_KEYVALUE,  # key=value形式として処理
        "spot_cation": OPTION_TYPE_KEYVALUE,  # key=value形式として処理
    }

    # 後処理関数で型変換
    def to_int(x):
        """文字列または数値を整数に変換"""
        if isinstance(x, int):
            return x
        return int(x)

    def to_int_tuple(x):
        """タプルまたはリストの各要素を整数に変換"""
        if isinstance(x, (tuple, list)):
            return tuple(int(v) for v in x)
        return (int(x),)

    def to_bool(x):
        """文字列または数値をブール値に変換"""
        if isinstance(x, bool):
            return x
        if isinstance(x, str):
            return x.lower() in ("true", "1", "yes", "on")
        return bool(x)

    post_processors = {
        "seed": to_int,
        "depol_loop": to_int,
        "replication_factors": to_int_tuple,
        "assess_cages": to_bool,
        # debugはフラグ型なのでpost_processor不要
    }

    processed, unprocessed = parse_options_generic(
        options, option_specs, post_processors
    )
    # processed（型変換済み）とunprocessed（プラグインで処理される可能性）の両方を統合
    return {**processed, **unprocessed}


# ============================================================================
# ヘルプ文章の定義
# ============================================================================

HELP_DEBUG = "Enable debug mode"
HELP_SHIFT = (
    "Shift the unit cell by the specified fractional coordinates along a, b, and c axes. "
    "This shifts all atomic positions by the given vector in fractional coordinates."
)
HELP_DEPOL_LOOP = (
    "Number of iterations for the depolarization optimization loop. "
    "Larger values may improve the quality of the hydrogen bond network. "
    "Default is 1000."
)
HELP_REPLICATION_MATRIX = (
    "Replication matrix as 9 integers (3x3 matrix). "
    "This matrix defines how the unit cell is replicated to create the supercell. "
    "The first row (p, q, r) specifies that the new a' axis direction is represented as pa+qb+rc "
    "using the original unit cell's basis vectors (a, b, c). "
    "Similarly, the second row (s, t, u) specifies that the b' axis direction is sa+tb+uc, "
    "and the third row defines the c' axis. "
    "For example, --replication_matrix 0 1 0  1 0 0  0 0 1 swaps the a and b axes of the unit cell. "
    "Another example, --replication_matrix 1 1 0  1 -1 0  0 0 1 transforms the unit cell such that "
    "a'=a+b and b'=a-b. "
    "If not specified, replication_factors is used instead."
)
HELP_REPLICATION_FACTORS = (
    "Repeat the unit cell along a, b, and c axes. "
    "For example, --rep 2 2 2 creates a 2x2x2 supercell. "
    "This is a convenient shortcut for diagonal replication matrices."
)
HELP_ASSESS_CAGES = (
    "Assess and report cage positions and types in the ice structure. "
    "This is useful for understanding the clathrate structure and "
    "determining where guest molecules can be placed."
)
HELP_GUEST = (
    "Guest molecule descriptors. Format: CAGE_TYPE=MOLECULE[*OCCUPANCY][+...]. "
    "Examples: A12=me (methane in A12 cages), A14=co2*0.5+et*0.3 "
    "(CO2 at 50% occupancy and ethane at 30% in A14 cages). "
    "Multiple guests can be specified with multiple -g options."
)
HELP_ANION = (
    "Specify a monatomic anion that replaces a water molecule. "
    "Format: LABEL=ION_NAME, where LABEL is the index of the water molecule in the unit cell. "
    "Examples: -a 1=Cl (replace water at position 1 with Cl-), -a 35=Br (replace at position 35 with Br-). "
    "Multiple anions can be specified with multiple -a options."
    "Note that they are replicated in the same way as water molecules when constructing the replicated cell."
)
HELP_CATION = (
    "Specify a monatomic cation that replaces a water molecule. "
    "Format: LABEL=ION_NAME, where LABEL is the index of the water molecule in the unit cell. "
    "Examples: -c 1=Na (replace water at position 1 with Na+), -c 35=K (replace at position 35 with K+). "
    "Multiple cations can be specified with multiple -c options."
    "Note that they are replicated in the same way as water molecules when constructing the replicated cell."
)
HELP_DENSITY = (
    "Target density of the ice in g/cm³. "
    "The unit cell will be scaled to achieve this density. "
    "If not specified, the original density of the lattice structure is used."
)
HELP_SEED = (
    "Random seed for guest molecule placement and other stochastic processes. "
    "Using the same seed will produce reproducible results."
)
HELP_SPOT_GUEST = (
    "Specify guest in the specific cage. "
    "Format: CAGE_INDEX=MOLECULE, where CAGE_INDEX is the index of the cage and MOLECULE is the guest molecule. "
    "Examples: -G 13=me (place methane in cage 13), -G 32=co2*0.5+et*0.3 (place CO2 at 50% occupancy and ethane at 30% in cage 32). "
    "Multiple spot guests can be specified with multiple -G options."
    "When the option argument is entered as ? (question mark), detailed information such as the position, type, and number of water molecules constituting the cage is displayed."
)
HELP_SPOT_ANION = (
    "Specify anion replacing the specified water molecule. "
    "Format: WATER_INDEX=ION_NAME, where WATER_INDEX is the index of the water molecule and ION_NAME is the anion name. "
    "Examples: -a 13=Cl (place Cl- in cage 13), -a 32=Br (place Br- in cage 32). "
    "Multiple spot anions can be specified with multiple -a options."
)
HELP_SPOT_CATION = (
    "Specify cation replacing the specified water molecule. "
    "Format: WATER_INDEX=ION_NAME, where WATER_INDEX is the index of the water molecule and ION_NAME is the cation name. "
    "Examples: -c 13=Na (place Na+ in cage 13), -c 32=K (place K+ in cage 32). "
    "Multiple spot cations can be specified with multiple -c options."
)
HELP_CONFIG = (
    "Path to a YAML configuration file. "
    "Settings from the config file will be overridden by command-line arguments. "
    "See documentation for the config file format."
)


def print_help():
    """ヘルプメッセージを表示"""
    print("Usage: genice3 [OPTIONS] UNITCELL")
    print()
    print("Options:")
    print(f"  -h, --help                Show this message and exit.")
    print(f"  --version, -V             {HELP_VERSION}")
    print(f"  -D, --debug               {HELP_DEBUG}")
    print(f"  --depol_loop INTEGER      {HELP_DEPOL_LOOP}")
    print(f"  --replication_matrix INT INT INT INT INT INT INT INT INT")
    print(f"                           {HELP_REPLICATION_MATRIX}")
    print(f"  --rep, --replication_factors INT INT INT")
    print(f"                           {HELP_REPLICATION_FACTORS}")
    print(f"  -s, --seed INTEGER        {HELP_SEED}")
    print(
        f"  -e, --exporter TEXT       Exporter plugin name (e.g., 'gromacs' or 'gromacs[options]')"
    )
    print(f"  -A, --assess_cages        {HELP_ASSESS_CAGES}")
    print(f"  -a, --spot_anion TEXT     {HELP_SPOT_ANION}")
    print(f"  -c, --spot_cation TEXT    {HELP_SPOT_CATION}")
    print(f"  -C, --config PATH         {HELP_CONFIG}")
    print()
    print("Arguments:")
    print("  UNITCELL                  Unitcell plugin name (required)")


def main() -> None:
    """メイン関数"""
    # --helpと--versionを先に処理
    if "--help" in sys.argv or "-h" in sys.argv:
        print_help()
        sys.exit(0)
    if "--version" in sys.argv or "-V" in sys.argv:
        print(f"genice3 {get_version()}")
        sys.exit(0)

    # ロギングを初期化（デフォルトはINFOレベル）
    basicConfig(level=INFO)
    logger = getLogger()

    # PoolBasedParserを使ってコマンドライン引数をパース
    parser = PoolBasedParser()
    try:
        parser.parse_args(sys.argv[1:])
    except Exception as e:
        logger.error(f"パースエラー: {e}")
        import traceback

        traceback.print_exc()
        sys.exit(1)

    result = parser.get_result()

    # バリデーション
    is_valid, errors = parser.validate()
    if not is_valid:
        for error in errors:
            logger.error(error)
        sys.exit(1)

    # 基底レベルのオプションを取得
    base_options = result["base_options"]
    # debugオプションが指定されていればloggerのlevelをDEBUGにする
    if base_options.get("debug"):
        logger.setLevel(DEBUG)
        for handler in logger.handlers:
            handler.setLevel(DEBUG)

    seed = base_options.get("seed", 1)
    depol_loop = base_options.get("depol_loop", 1000)
    assess_cages = base_options.get("assess_cages", False)
    spot_anion_dict = base_options.get("spot_anion", {}) or {}
    spot_cation_dict = base_options.get("spot_cation", {}) or {}

    # replication_matrixとreplication_factorsの処理
    replication_matrix = base_options.get("replication_matrix")
    replication_factors = base_options.get("replication_factors", (1, 1, 1))
    if replication_matrix is None:
        # replication_factorsがリストの場合はタプルに変換
        if isinstance(replication_factors, list):
            replication_factors = tuple(replication_factors)
        replication_matrix = np.diag(replication_factors)
    else:
        replication_matrix = np.array(replication_matrix)
    logger.debug(
        f"Settings: {seed=} {depol_loop=} {assess_cages=} {spot_anion_dict=} {spot_cation_dict=} {replication_matrix=}"
    )

    # indexを数字に変換する。
    spot_anion_dict = ion_processor(spot_anion_dict)
    spot_cation_dict = ion_processor(spot_cation_dict)

    # unitcellプラグインの処理結果
    unitcell_name = result["unitcell"]["name"]
    unitcell_processed = result["unitcell"]["processed"]

    # exporterプラグインの処理結果
    exporter_name = result["exporter"]["name"] or "gromacs"  # デフォルトはgromacs
    exporter_processed = result["exporter"]["processed"]

    # GenIce3インスタンスを作成
    genice = GenIce3(
        depol_loop=depol_loop,
        replication_matrix=replication_matrix,
        seed=seed,
        spot_anions=spot_anion_dict,
        spot_cations=spot_cation_dict,
    )

    # unitcellプラグインを設定
    genice.unitcell = safe_import("unitcell", unitcell_name).UnitCell(
        **unitcell_processed
    )

    # unitcellを設定したので、cage調査ができる。
    if assess_cages:
        survey_result = genice.cage_survey
        print(survey_result)
        sys.exit(0)

    # コマンドライン全体を取得
    command_line = " ".join(sys.argv)
    exporter_processed["command_line"] = command_line

    # exporterプラグインを実行
    safe_import("exporter", exporter_name).dump(
        genice, sys.stdout, **exporter_processed
    )


if __name__ == "__main__":
    main()
