import click
from logging import getLogger, basicConfig, DEBUG, INFO
from typing import Optional, Tuple, Dict, Any
import numpy as np
import sys
from typing import Dict, List
import json
from importlib.metadata import version, PackageNotFoundError

from genice3.plugin import parse_plugin_options, safe_import
from genice3.genice import GenIce3
from genice3.config import load_config, parse_config, merge_config_with_args
from genice3.unitcell import ion_processor

HELP_VERSION = "Show the version and exit."


def get_version():
    """パッケージのバージョンを取得する。見つからない場合はデフォルト値を返す。"""
    try:
        return version("genice3")
    except PackageNotFoundError:
        return "3.*.*.*"


# ============================================================================
# Clickヘルプ文章の定義
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


# def _ion_parser(ion_options: List[str]) -> Dict[int, str]:
#     """イオンオプションをパースする。常にDict[int, str]を返す。"""
#     from genice3.molecule.one import Molecule

#     result: Dict[int, Molecule] = {}
#     for ion_option in ion_options:
#         label, ion_name = ion_option.split("=")
#         result[int(label)] = Molecule(name=ion_name, label=ion_name)
#     return result


def _merge_hierarchical_dict(
    base: Dict[str, Any], override: Dict[str, Any]
) -> Dict[str, Any]:
    """
    階層構造の辞書をマージする（overrideが優先）

    Args:
        base: ベースとなる辞書（設定ファイルから）
        override: 上書きする辞書（コマンドライン引数から）

    Returns:
        マージされた辞書

    Examples:
        >>> base = {"guest": {"A12": "me"}, "shift": [0.1, 0.1, 0.1]}
        >>> override = {"guest": {"A14": "et"}}
        >>> _merge_hierarchical_dict(base, override)
        {"guest": {"A12": "me", "A14": "et"}, "shift": [0.1, 0.1, 0.1]}
    """
    result = base.copy()
    for key, value in override.items():
        if key in result and isinstance(result[key], dict) and isinstance(value, dict):
            # 両方とも辞書の場合は再帰的にマージ
            result[key] = _merge_hierarchical_dict(result[key], value)
        else:
            # それ以外の場合は上書き
            result[key] = value
    return result


# プラグインのオプションと設定を統合する。
def _plugin_settings(plugin_type: str, cmdline: str, config: dict, default=None):
    # pluginオプションを統合
    # configをcmdlineが上書きする
    logger = getLogger()

    name = ""
    options = None

    if config and isinstance(config, dict) and "name" in config:
        # 設定ファイルから取得
        name = config["name"]
        options = config
        logger.debug(f"config: {name} {options}")

    if cmdline:
        name, options = parse_plugin_options(cmdline)
        logger.debug(f"cmdline: {name} {options}")
    elif name == "":
        if default:
            name, options = parse_plugin_options(default)
            logger.debug(f"default: {name} {options}")
        else:
            raise ValueError(f"Invalid type for {plugin_type}: {type(cmdline)}")

    logger.debug(f"finally: {name} {options}")
    return name, options


# clickを用い、-dオプションでデバッグレベルを指定できるようにする。
@click.command()
@click.help_option()
@click.version_option(version=get_version(), prog_name="genice3", help=HELP_VERSION)
@click.argument("unitcell", type=str, required=False)
@click.option("--debug", "-D", is_flag=True, help=HELP_DEBUG)
@click.option("--depol_loop", type=int, default=None, help=HELP_DEPOL_LOOP)
@click.option(
    "--replication_matrix",
    type=click.Tuple([int, int, int, int, int, int, int, int, int]),
    default=None,
    help=HELP_REPLICATION_MATRIX,
)
@click.option(
    "--replication_factors",
    "--rep",
    type=click.Tuple([int, int, int]),
    default=None,
    help=HELP_REPLICATION_FACTORS,
)
@click.option("--seed", "-s", type=int, default=None, help=HELP_SEED)
@click.option(
    "--exporter",
    "-e",
    default=None,
    help="Exporter plugin name (e.g., 'gromacs' or 'gromacs[options]')",
)
# assess_cagesをGenIce3のオプションに戻す。
@click.option(
    "--assess_cages", "-A", is_flag=True, default=False, help=HELP_ASSESS_CAGES
)
@click.option("--spot_anion", "-a", multiple=True, default=None, help=HELP_SPOT_ANION)
@click.option("--spot_cation", "-c", multiple=True, default=None, help=HELP_SPOT_CATION)
@click.option(
    "--config", "-C", type=click.Path(exists=True), default=None, help=HELP_CONFIG
)
def main(
    unitcell: Optional[str],
    debug: bool,
    depol_loop: int,
    replication_matrix: Optional[Tuple[int, int, int, int, int, int, int, int, int]],
    replication_factors: Tuple[int, int, int],
    seed: int,
    exporter: str,
    assess_cages: Optional[bool],
    spot_anion: List[str],
    spot_cation: List[str],
    config: Optional[str],
) -> None:
    basicConfig(level=DEBUG if debug else INFO)
    logger = getLogger()

    # 設定ファイルを読み込む（存在する場合）
    config_options = {}
    if config:
        try:
            config_data = load_config(config)
            logger.info(f"設定ファイルから読み込んだデータ: {config_data}")
            config_options = parse_config(config_data)
            logger.info(f"設定ファイルから読み込んだオプション: {config_options}")
        except Exception as e:
            logger.warning(f"設定ファイルの読み込みに失敗しました: {e}")
            logger.warning("設定ファイルを無視して続行します。")

    # コマンドライン引数を辞書に変換
    args_dict = {
        "seed": seed,
        "depol_loop": depol_loop,
        "replication_matrix": replication_matrix,
        "replication_factors": replication_factors,
        "assess_cages": assess_cages,
        "debug": debug,
        "spot_anion": dict([keyvalue.split("=") for keyvalue in spot_anion]),
        "spot_cation": dict([keyvalue.split("=") for keyvalue in spot_cation]),
        "exporter": exporter,
        "unitcell": unitcell,
    }

    # 設定ファイルからunitcellとexporterセクションを直接取得（プラグインにそのまま渡す）
    unitcell_config = config_options.get("unitcell")
    exporter_config = config_options.get("exporter")

    if spot_anion:
        spot_anion = dict([keyvalue.split("=") for keyvalue in spot_anion])
    else:
        spot_anion = None
    if spot_cation:
        spot_cation = dict([keyvalue.split("=") for keyvalue in spot_cation])
    else:
        spot_cation = None

    # 設定ファイルとコマンドライン引数を統合（コマンドライン引数が優先）
    # unitcellとexporterは除外（プラグイン側で処理するため）
    args_dict_without_plugins = {
        "seed": seed,
        "depol_loop": depol_loop,
        "replication_matrix": replication_matrix,
        "replication_factors": replication_factors,
        "assess_cages": assess_cages,
        "debug": debug,
        "spot_anion": spot_anion,
        "spot_cation": spot_cation,
    }
    merged_options = merge_config_with_args(config_options, args_dict_without_plugins)

    # 統合されたオプションから値を取得
    seed = merged_options.get("seed", 1)
    depol_loop = merged_options.get("depol_loop", 1000)
    assess_cages = merged_options.get("assess_cages", False)
    debug = merged_options.get("debug", False)
    spot_anion_dict = merged_options.get("spot_anion", {}) or {}
    spot_cation_dict = merged_options.get("spot_cation", {}) or {}
    logger.debug(
        f"Raw settings for main(): {seed=} {depol_loop=} {assess_cages=} {debug=} {spot_anion_dict=} {spot_cation_dict=} {replication_matrix=} {replication_factors=}"
    )

    # indexを数字に変換する。
    spot_anion_dict = ion_processor(spot_anion_dict)
    spot_cation_dict = ion_processor(spot_cation_dict)

    # replication_matrixとreplication_factorsの処理
    replication_matrix = merged_options.get("replication_matrix")
    replication_factors = merged_options.get("replication_factors", (1, 1, 1))
    if replication_matrix is None:
        replication_matrix = np.diag(replication_factors)
    else:
        replication_matrix = np.array(replication_matrix)
    logger.debug(
        f"Cooked settings for main(): {seed=} {depol_loop=} {assess_cages=} {debug=} {spot_anion_dict=} {spot_cation_dict=} {replication_matrix=}"
    )

    exporter_name, exporter_options = _plugin_settings(
        "exporter", exporter, exporter_config, default="gromacs"
    )
    unitcell_name, unitcell_options = _plugin_settings(
        "unitcell", unitcell, unitcell_config
    )

    logger.debug("Debug mode enabled")
    genice = GenIce3(
        depol_loop=depol_loop,
        replication_matrix=replication_matrix,
        seed=seed,
        spot_anions=spot_anion_dict,
        spot_cations=spot_cation_dict,
    )
    # genice.unitcell = Ice1h(shift=shift, assess_cages=assess_cages)
    genice.unitcell = safe_import("unitcell", unitcell_name).UnitCell(
        **unitcell_options
    )

    # unitcellを設定したので、cage調査ができる。
    if assess_cages:
        survey_result = genice.cage_survey
        print(survey_result)
        sys.exit(0)

    # コマンドライン全体を取得
    command_line = " ".join(sys.argv)
    exporter_options["command_line"] = command_line
    safe_import("exporter", exporter_name).dump(genice, sys.stdout, **exporter_options)


if __name__ == "__main__":
    main()
