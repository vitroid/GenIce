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


def compose_long_line_from_list(items: list, N=70):
    while True:
        line = ""
        while len(line) < N and items:
            line += " " + items.pop(0)
        yield line + "\\\n"
        if not items:
            yield line + "\\\n"
            break


basicConfig(level=INFO)
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


def make_test(ice: str, tests: dict, formatters: dict):
    for formatter_prefix, formatter_path in formatters.items():
        for i, test in enumerate(tests):
            product = f"{ice}_{i}.{formatter_prefix}"

            genice_options = ""
            module_options = ""

            logger.debug(f"{ice} {test}")
            if "options" in test:
                genice_options = test["options"]
            if "args" in test:
                module_options = options_parser(test["args"])

            genice_options += " " + " ".join(random.sample(additional_options, 3))
            genice_options += f" -f {formatter_prefix}"
            if module_options != "":
                module_options = "[" + module_options + "]"
            target = f"{ice}{module_options}"
            logger.debug(f"Target: {target}")
            # if testmode:
            rule = f"{product}.diff: {product} ../../genice2/lattices/{ice}.py {formatter_path}\n"
            rule += f"\t$(GENICE) {target} {genice_options} | diff - $<\n\ttouch $@\n"
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
            formatter_paths[formatter_prefix] = filepath
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

    for product, rule in make_test(ice, tests[ice], formatters):
        products_prepare.append(product)
        products_test.append(product)
        rules += rule

print("GENICE=../../genice.x")
# print("GENICE=genice2")
targets = compose_long_line_from_list(products_prepare)
print("TARGETS=", *targets)
print("prepare: $(TARGETS)")
print(
    "test: $(foreach file,$(TARGETS),$(if $(shell test -s $(file) && echo 1),$(file).diff,))"
)
print(rules)
