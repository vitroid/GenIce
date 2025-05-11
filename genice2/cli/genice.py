import sys
import argparse as ap
from collections import defaultdict
import logging
from typing import Optional

from genice2.plugin import safe_import, descriptions

# from genice2 import __version__
from genice2.cli import SmartFormatter, logger_setup, help_water, help_format
from genice2.genice import GenIce, GenIceConfig
from genice2.valueparser import plugin_option_parser, parse_guest
from genice2.decorators import timeit, banner
import pickle
import numpy as np
from importlib.metadata import version, PackageNotFoundError
import click
from genice2.genice import GenIce, GenIceConfig

from genice2.lattices import Lattice
from genice2.molecules import Molecule
from genice2.formats import Format

# from genice2.water import water_models
# from genice2.lattices import lattices
# from genice2.decorators import SmartFormatter

# 遅延評価。descriptions()関数は重いので、必要なければ呼びたくない。

# ヘルプメッセージの定数定義
HELP_VERSION = "Show the version and exit."
HELP_REP = "Repeat the unit cell along a, b, and c axes. [1,1,1]"
HELP_RESHAPE = (
    "Convert the unit cell shape by specifying the new (a,b,c) set from the original (a,b,c) of the unit cell. "
    "The combination of (a,b,c) is specified by nine integers. "
    "For example, '--reshape 3,0,0,0,2,0,0,0,1' specifies that the new cell vectors are (3a, 2b, c), "
    "which is equivalent to '--rep 3 2 1'."
)
HELP_SHIFT = "Shift the unit cell along a, b, and c axes. (0.5==half cell) [0,0,0]"
HELP_DENS = "Specify the ice density in g/cm3 (Guests are not included.)"
HELP_NOISE = (
    "Add a Gauss noise with given width (SD) to the molecular positions of water. "
    "The value 1 corresponds to 1 percent of the molecular diameter of water."
)
HELP_SEED = "Random seed [1000]"
HELP_FORMAT = "Output format"
HELP_WATER = "Water model"
HELP_GUEST = (
    "Specify guest(s) in the cage type. (D=empty, T=co2*0.5+me*0.3, etc.)\n\n"
    + descriptions("molecule", water=False, width=55)
)
HELP_SPOT_GUEST = "Specify guest in the specific cage. (13=me, 32=co2, etc.)"
HELP_GROUP = "Specify the group. (-H 13=bu-:0, etc.)"
HELP_ANION = (
    "Specify a monatomic anion that replaces a water molecule. (3=Cl, 39=F, etc.)"
)
HELP_CATION = (
    "Specify a monatomic cation that replaces a water molecule. (3=Na, 39=NH4, etc.)"
)
HELP_DEPOL = 'Depolarization. (strict, optimal, or none) ["strict"]'
HELP_ASIS = "Assumes all given HB pairs to be fixed. No shuffle and no depolarization."
HELP_DEBUG = "Output debugging info."
HELP_QUIET = "Do not output progress messages."
HELP_ASSESS_CAGES = (
    "Assess the locations of cages based on the HB network topology. "
    "Note: it may fail when the unit cell is too small."
)
HELP_OUTPUT = "Output file name"
# HELP_TYPE = (
#     "R|Crystal type (1c, 1h, etc. See https://github.com/vitroid/GenIce for available ice structures.)\n\n"
#     "If you want to analyze your own structures, please try analice tool.\n\n"
#     + descriptions("lattice", width=55)
# )


def help_type():
    return (
        "R|Crystal type (1c, 1h, etc. See https://github.com/vitroid/GenIce for available ice structures.)\n\nIf you want to analyze your own structures, please try analice tool.\n\n"
        + descriptions("lattice", width=55)
    )


def help_guest():
    return (
        "R|Specify guest(s) in the cage type. (D=empty, T=co2*0.5+me*0.3, etc.)\n\n"
        + descriptions("molecule", water=False, width=55)
    )


def get_version():
    """パッケージのバージョンを取得する。見つからない場合はデフォルト値を返す。"""
    try:
        return version("genice2")
    except PackageNotFoundError:
        return "2.*.*.*"


def setup_logging(debug: bool, quiet: bool) -> logging.Logger:
    """ログレベルの設定を行う"""
    logger = logging.getLogger()
    if debug:
        logging.basicConfig(level=logging.DEBUG)
    elif quiet:
        logging.basicConfig(level=logging.WARNING)
    else:
        logging.basicConfig(level=logging.INFO)
    return logger


def setup_water_model(water_spec: str, logger: logging.Logger) -> Molecule:
    """水分子モデルの設定を行う"""
    water_type, water_options = plugin_option_parser(water_spec)
    logger.debug(f"Water type: {water_type}")
    return safe_import("molecule", water_type).Molecule(**water_options)


def setup_lattice(ice_type: str, logger: logging.Logger) -> Lattice:
    """格子の設定を行う"""
    lattice_type, lattice_options = plugin_option_parser(ice_type)
    logger.debug(f"Lattice: {lattice_type}")
    assert lattice_type is not None
    return safe_import("lattice", lattice_type).Lattice(**lattice_options)


def setup_formatter(format_spec: str, logger: logging.Logger) -> Format:
    """フォーマッターの設定を行う"""
    file_format, format_options = plugin_option_parser(format_spec)
    logger.debug(f"Output file format: {file_format}")
    formatter_module = safe_import("format", file_format)
    return formatter_module.Format(**format_options)


def create_config(kwargs: dict, water: Molecule) -> GenIceConfig:
    """設定オブジェクトの作成"""
    config = GenIceConfig(
        water=water,
        density=kwargs["dens"],
        rep=tuple(kwargs["rep"]),
        reshape=parse_reshape(kwargs["reshape"]),
        shift=tuple(kwargs["shift"]),
        noise=kwargs["add_noise"],
        depol=kwargs["depol"],
        asis=kwargs["asis"],
        assess_cages=kwargs["assess_cages"],
    )

    # ゲスト分子の設定
    if kwargs["guest"]:
        config.guests = parse_guests(kwargs["guest"])
    if kwargs["spot_guest"]:
        config.spot_guests = parse_spot_guests(kwargs["spot_guest"])
    if kwargs["group"]:
        config.spot_groups = parse_groups(kwargs["group"])
    if kwargs["anion"]:
        config.anions = parse_ions(kwargs["anion"])
    if kwargs["cation"]:
        config.cations = parse_ions(kwargs["cation"])

    return config


def output_result(ice: GenIce, formatter: Format, output_file: Optional[str] = None):
    """結果の出力"""
    if output_file:
        with open(output_file, "w") as f:
            formatter.dump(ice, f)
    else:
        formatter.dump(ice, sys.stdout)


@click.command()
@click.version_option(version=get_version(), prog_name="genice2", help=HELP_VERSION)
@click.option("--rep", "-r", nargs=3, type=int, default=[1, 1, 1], help=HELP_REP)
@click.option("--reshape", "-R", type=str, default="", help=HELP_RESHAPE)
@click.option(
    "--shift", "-S", nargs=3, type=float, default=[0.0, 0.0, 0.0], help=HELP_SHIFT
)
@click.option("--dens", "-d", type=float, default=-1, help=HELP_DENS)
@click.option(
    "--add_noise", type=float, default=0.0, metavar="percent", help=HELP_NOISE
)
@click.option("--seed", "-s", type=int, default=1000, help=HELP_SEED)
@click.option("--format", "-f", default="gromacs", metavar="name", help=HELP_FORMAT)
@click.option("--water", "-w", default="tip3p", metavar="model", help=HELP_WATER)
@click.option("--guest", "-g", multiple=True, metavar="14=me", help=HELP_GUEST)
@click.option(
    "--spot_guest", "-G", multiple=True, metavar="13=me", help=HELP_SPOT_GUEST
)
@click.option("--group", "-H", multiple=True, metavar="13=bu-:0", help=HELP_GROUP)
@click.option("--anion", "-a", multiple=True, metavar="3=Cl", help=HELP_ANION)
@click.option("--cation", "-c", multiple=True, metavar="3=Na", help=HELP_CATION)
@click.option("--depol", default="strict", help=HELP_DEPOL)
@click.option("--asis", is_flag=True, default=False, help=HELP_ASIS)
@click.option("--debug", "-D", is_flag=True, help=HELP_DEBUG)
@click.option("--quiet", "-q", is_flag=True, help=HELP_QUIET)
@click.option("--assess_cages", "-A", is_flag=True, help=HELP_ASSESS_CAGES)
@click.option("--output", "-o", type=str, help=HELP_OUTPUT)
@click.argument("ice_type", type=str)
def main_click(ice_type, **kwargs):
    """GenIce is a swiss army knife to generate hydrogen-disordered ice structures."""
    # ログレベルの設定
    logger = setup_logging(kwargs["debug"], kwargs["quiet"])

    # ランダムシードの設定
    if kwargs["seed"] is not None:
        np.random.seed(kwargs["seed"])

    # 水分子モデルの設定
    water = setup_water_model(kwargs["water"], logger)

    # 格子の設定
    lattice = setup_lattice(ice_type, logger)

    # フォーマッターの設定
    formatter = setup_formatter(kwargs["format"], logger)

    # 設定オブジェクトの作成
    config = create_config(kwargs, water)

    # 氷の生成
    ice = GenIce(lattice, config)

    # 結果の出力
    output_result(ice, formatter, kwargs.get("output"))


def parse_reshape(reshape_str):
    if not reshape_str:
        return np.eye(3, dtype=int)
    try:
        return np.array([int(x) for x in reshape_str.split(",")]).reshape(3, 3)
    except:
        raise click.BadParameter("reshape must be 9 integers separated by commas")


def parse_guests(guests):
    result = {}
    for guest in guests:
        try:
            cage, molecule = guest.split("=")
            result[int(cage)] = molecule
        except:
            raise click.BadParameter('guest must be in the format "cage=molecule"')
    return result


def parse_spot_guests(guests):
    result = {}
    for guest in guests:
        try:
            cage, molecule = guest.split("=")
            result[int(cage)] = molecule
        except:
            raise click.BadParameter('spot_guest must be in the format "cage=molecule"')
    return result


def parse_groups(groups):
    result = {}
    for group in groups:
        try:
            cage, group_info = group.split("=")
            result[int(cage)] = group_info
        except:
            raise click.BadParameter('group must be in the format "cage=group_info"')
    return result


def parse_ions(ions):
    result = {}
    for ion in ions:
        try:
            pos, element = ion.split("=")
            result[int(pos)] = element
        except:
            raise click.BadParameter('ion must be in the format "position=element"')
    return result


if __name__ == "__main__":
    main_click()
