# test of tests
import sys
import random
import glob
import os

sys.path.append("../..")
from genice2.plugin import scan
from logging import getLogger, basicConfig, DEBUG, INFO


def options_parser(options):
    if type(options) == dict:
        options = ":".join([f"{x}={y}" for x, y in options.items()])
    return options


basicConfig(level=DEBUG)
logger = getLogger()

result = scan("lattice")
# preinstalled (system) plugins only
ices = result["system"]
tests = result["tests"]
# print(tests)
products_prepare = []
products_test = []
rules = []

additional_options = [
    "--reshape 1,1,1,1,-1,0,1,1,-2",
    # "--assess_cages",
    "--shift 0.5 0.4 0.3",
    "--add_noise 0.01",
    "--water 4site",
    "--cation 0=Na --anion 2=Cl",
    # "--asis",
]
formatter_paths = {}
for filepath in glob.glob("../../genice2/formats/*.py"):
    if not os.path.islink(filepath):
        formatter_prefix = os.path.basename(filepath).split(".")[0]
        if formatter_prefix == "raw":
            continue
        formatter_paths[f" -f {formatter_prefix}"] = filepath

for ice in ices:
    if ice in ("iceR",):
        continue
    if ice not in tests:
        tests[ice] = ("",)

    for i, test in enumerate(tests[ice]):
        genice_options = ""
        module_options = ""

        if type(test) is dict:
            logger.debug(f"{ice} {test}")
            if "options" in test:
                genice_options = test["options"]
            if "args" in test:
                module_options = options_parser(test["args"])
        # other available options:
        # - noise: float
        # - assess_cages: bool
        # - depol: str
        # - shift: tuple
        # - rep: tuple
        # - reshape: np.ndarray
        # - signature: str
        # - density: float
        # - cations: dict
        # - anions: dict
        # - spot_guests: dict
        # - spot_groups: dict
        # - asis: bool

        genice_options += " " + " ".join(random.sample(additional_options, 3))
        format_option = random.choice(list(formatter_paths.keys()))
        genice_options += format_option
        product = f"{ice}_{i}.output"
        if module_options != "":
            module_options = "[" + module_options + "]"
        target = f"{ice}{module_options}"
        logger.debug(f"Target: {target}")
        # if testmode:
        rules.append(
            f"{product}.diff: {product} ../../genice2/lattices/{ice}.py {formatter_paths[format_option]}\n"
        )
        rules.append(
            f"\t$(GENICE) {target} {genice_options} | diff - $< && touch $@\n"
        )  # > $@.diff\n")
        products_test.append(f"{product}.diff")
        # else:
        rules.append(
            f"{product}: ../../genice2/lattices/{ice}.py  {formatter_paths[format_option]}\n"
        )
        rules.append(f"\t$(GENICE) {target} {genice_options} > $@\n")
        products_prepare.append(product)

print("GENICE=../../genice.x")
# print("GENICE=genice2")
print("TARGETS=", *products_prepare)
print("prepare: $(TARGETS)")
print("test: $(patsubst %, %.diff, $(TARGETS))")
print("".join(rules))
