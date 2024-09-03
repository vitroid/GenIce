# test of tests
import sys

sys.path.append("/Users/matto/Unison/github/GenIce2")
from genice2.plugin import scan
from logging import getLogger, basicConfig, DEBUG, INFO


def options_parser(options):
    if type(options) == dict:
        options = ":".join([f"{x}={y}" for x, y in options.items()])
    return options


basicConfig(level=DEBUG)
logger = getLogger()

testmode = len(sys.argv) > 1 and sys.argv[1] == "test"

result = scan("lattice")
# preinstalled (system) plugins only
ices = result["system"]
tests = result["tests"]
# print(tests)
products = []
rules = []
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
        product = f"{ice}_{i}.gro"
        if module_options != "":
            module_options = "[" + module_options + "]"
        target = f"{ice}{module_options}"
        logger.debug(f"Target: {target}")
        if testmode:
            rules.append(f"{product}.diff: {product} ../../genice2/lattices/{ice}.py\n")
            rules.append(
                f"\t../../genice.x {target} {genice_options} | diff - $< \n"
            )  # > $@.diff\n")
            products.append(f"{product}.diff")
        else:
            rules.append(f"{product}: ../../genice2/lattices/{ice}.py\n")
            rules.append(f"\t../../genice.x {target} {genice_options} > $@\n")
            products.append(product)

print("all:", *products)
print("".join(rules))
