import click
from logging import getLogger, basicConfig, DEBUG, INFO
from typing import Optional, Tuple
import numpy as np
import sys
from typing import Dict, List
import json


from genice2.cli.genice import get_version, HELP_VERSION

from genice3.plugin import parse_plugin_specification, safe_import
from genice3.genice import GenIce3


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


def _ion_parser(ion_options: List[str]) -> Dict[int, str]:
    """イオンオプションをパースする。常にDict[int, str]を返す。"""
    result: Dict[int, str] = {}
    for ion_option in ion_options:
        label, ion_name = ion_option.split("=")
        result[int(label)] = ion_name
    return result


# clickを用い、-dオプションでデバッグレベルを指定できるようにする。
@click.command()
@click.help_option()
@click.version_option(version=get_version(), prog_name="genice3", help=HELP_VERSION)
@click.argument("unitcell", type=str)
@click.option("--debug", "-D", is_flag=True, help=HELP_DEBUG)
@click.option("--depol_loop", type=int, default=1000, help=HELP_DEPOL_LOOP)
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
    default=(1, 1, 1),
    help=HELP_REPLICATION_FACTORS,
)
@click.option("--seed", "-s", type=int, default=1, help=HELP_SEED)
@click.option("--exporter", "-e", default="gromacs")
# assess_cagesをGenIce3のオプションに戻す。
@click.option("--assess_cages", "-A", is_flag=True, help=HELP_ASSESS_CAGES)
def main(
    unitcell: str,
    debug: bool,
    depol_loop: int,
    replication_matrix: Optional[Tuple[int, int, int, int, int, int, int, int, int]],
    replication_factors: Tuple[int, int, int],
    seed: int,
    exporter: str,
    assess_cages: bool,
) -> None:
    basicConfig(level=DEBUG if debug else INFO)
    logger = getLogger()
    if replication_matrix is None:
        replication_matrix = np.diag(replication_factors)
    else:
        replication_matrix = np.array(replication_matrix)

    # exporterオプションをパース：プラグイン名とパース済みオプション辞書を取得
    exporter_name, exporter_options = parse_plugin_specification(exporter)
    exporter_module = safe_import("exporter", exporter_name)

    # unitcellオプションをパース
    unitcell_name, unitcell_options = parse_plugin_specification(unitcell)
    unitcell_module = safe_import("unitcell", unitcell_name)

    # # cage?フラグがある場合、GenIce3の作成に必要な情報をunitcell_optionsに追加
    # # （UnitCell.__init__内でデフォルト値が使われるが、CLI側の値を使うため）
    # if "cage?" in unitcell_options:
    #     unitcell_options["replication_matrix"] = replication_matrix
    #     unitcell_options["depol_loop"] = depol_loop
    #     unitcell_options["seed"] = seed

    logger.debug("Debug mode enabled")
    genice = GenIce3(
        depol_loop=depol_loop,
        replication_matrix=replication_matrix,
        seed=seed,
    )
    logger.info("Reactive properties:")
    logger.info(f"     All: {genice.list_all_reactive_properties().keys()}")
    logger.info(f"  Public: {genice.list_public_reactive_properties().keys()}")
    logger.info("Settabe reactive properties:")
    logger.info(f"     All: {genice.list_settable_reactive_properties().keys()}")
    logger.info(f"  Public: {genice.list_public_settable_reactive_properties().keys()}")
    # genice.unitcell = Ice1h(shift=shift, assess_cages=assess_cages)
    genice.unitcell = unitcell_module.UnitCell(**unitcell_options)

    # unitcellを設定したので、cage調査ができる。
    if assess_cages:
        survey_result = genice.cage_survey
        print(survey_result)
        sys.exit(0)
    # # stage 2
    # print(genice.graph)
    # # stage 4
    # print(genice.digraph)
    # # genice.unitcell = ice1h(shift=shift)
    # # stage 5
    # # print(genice.orientations)
    # # stage 6 and 7
    # print(genice.molecules(types=[MoleculeType.WATER]))
    # waters = genice.water_molecules(water_model=FourSiteWater())
    # guests = genice.guest_molecules(guests=guest_info, spot_guests=spot_guest_info)
    # ions = genice.substitutional_ions()
    # spot_dopants: これも要るだろうね。
    # spot_guest: その前に、replicate後のcageの場所を何らかの方法でユーザーに知らせないといけない。まああとまわしでいいだろう。
    # group

    # 試しにgroファイルにしてみるか。
    # with open("genice3.gro", "w") as f:
    #     f.write(to_gro(cellmat=genice.cell, waters=waters, guests=guests, ions=ions))

    # コマンドライン全体を取得
    command_line = " ".join(sys.argv)

    # パース済みのオプション辞書をプラグインに渡す
    # プラグイン側のデコレータは文字列をパースするが、辞書が渡された場合はそのまま通す
    # logger.info(f"{exporter_options=}")
    exporter_options["command_line"] = command_line
    exporter_module.dump(genice, sys.stdout, **exporter_options)


if __name__ == "__main__":
    main()
