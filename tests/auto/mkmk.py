# test of tests
import sys
import random
import glob
import os

sys.path.append("../..")
from genice2 import plugin
from logging import getLogger, basicConfig, DEBUG, INFO


def options_parser(options):
    if type(options) == dict:
        options = ":".join([f"{x}={y}" for x, y in options.items()])
    return options


basicConfig(level=DEBUG)
logger = getLogger()

# print(tests)


# other available options:
# - noise: float
# - assess_cages: bool
# - depol: str
## - shift: tuple
# - rep: tuple
## - reshape: np.ndarray
# - signature: str
# - density: float
# - cations: dict
# - anions: dict
# - spot_guests: dict
# - spot_groups: dict
# - asis: bool


def make_test(ice: str, tests: dict, format: str, formatter_path: str):
    if ice == "HS1":
        logger.info(
            f"Make tests of {ice} tests {tests} in {format} via {formatter_path}"
        )

    for i, test in enumerate(tests):
        product = f"{ice}_{i}.output"

        genice_options = ""
        module_options = ""

        logger.debug(f"{ice} {test}")
        if "options" in test:
            genice_options = test["options"]
        if "args" in test:
            module_options = options_parser(test["args"])

        genice_options += " " + " ".join(random.sample(additional_options, 3))
        genice_options += format
        if module_options != "":
            module_options = "[" + module_options + "]"
        target = f"{ice}{module_options}"
        logger.debug(f"Target: {target}")
        # if testmode:
        rule = f"{product}.diff: {product} ../../genice2/lattices/{ice}.py {formatter_path}\n"
        rule += f"\t$(GENICE) {target} {genice_options} | diff - $< && touch $@\n"
        rule += f"{product}: ../../genice2/lattices/{ice}.py  {formatter_path}\n"
        rule += f"\t$(GENICE) {target} {genice_options} > $@\n"
        yield product, rule


additional_options = [
    "--reshape 1,1,1,1,-1,0,1,1,-2",
    # "--assess_cages",
    "--shift 0.5 0.4 0.3",
    "--add_noise 0.01",
    "--water 4site",
    "--cation 0=Na --anion 2=Cl",
    # "--asis",
]


def formatter_list():
    formatter_paths = {}
    for filepath in glob.glob("../../genice2/formats/*.py"):
        if not os.path.islink(filepath):
            formatter_prefix = os.path.basename(filepath).split(".")[0]
            if formatter_prefix in ("raw", "null", "__init__"):
                continue
            formatter_paths[f" -f {formatter_prefix}"] = filepath
    return formatter_paths


result = plugin.scan("lattice")
# preinstalled (system) plugins only
ices = result["system"]
tests = result["tests"]
formatters = formatter_list()

# print(tests)
# sys.exit(0)

products_prepare = []
products_test = []
rules = ""


for ice in ices:
    if ice in ("iceR",):
        continue
    if ice not in tests:
        tests[ice] = ("",)
    format = random.choice(list(formatters.keys()))

    for product, rule in make_test(ice, tests[ice], format, formatters[format]):
        products_prepare.append(product)
        products_test.append(product)
        rules += rule

print("GENICE=../../genice.x")
# print("GENICE=genice2")
print("TARGETS=", *products_prepare)
print("prepare: $(TARGETS)")
print("test: $(patsubst %, %.diff, $(TARGETS))")
print(rules)
