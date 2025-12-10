from logging import getLogger
from typing import Dict, List
from genice3.genice import ConfigurationError, ShowUsageError, GuestSpec
from genice3.plugin import safe_import
from genice3.molecule import Molecule


def guest_processor(arg: dict):
    # 辞書のvaluesの特殊記法をparseする。
    logger = getLogger("guest_processor")
    logger.info(f"{arg=}")
    result = {}
    for cage, guest_specs in arg.items():
        logger.info(f"{cage=} {guest_specs=}")
        result[cage] = []
        total_occupancy = 0
        for guest_spec in guest_specs.split("+"):
            if "*" in guest_spec:
                occupancy, molecule = guest_spec.split("*")
                occupancy = float(occupancy)
            else:
                molecule = guest_spec
                occupancy = 1.0
            # moleculeはオプションを含んでいるかもしれないが、今は処理しない。
            molecule = safe_import("molecule", molecule).Molecule()
            result[cage].append(GuestSpec(molecule, occupancy))
            total_occupancy += occupancy
        if total_occupancy > 1.0:
            raise ConfigurationError(f"Total occupancy of {option} is greater than 1.0")
    logger.debug(f"{result=}")
    return result


def spot_guest_processor(arg: dict) -> Dict[int, Molecule]:
    # keyとvalueを変換するのみ
    result: Dict[int, Molecule] = {}
    for label, molecule in arg.items():
        if label == "?":
            raise ShowUsageError(
                "Use '?' to display cage information. "
                "This option triggers a survey of cage positions and types."
            )
        result[int(label)] = safe_import("molecule", molecule).Molecule()
    return result


def water_model_processor(arg: str) -> Molecule:
    # オプション文字列は今のところないのでこれで動く。
    return safe_import("molecule", arg).Molecule()
