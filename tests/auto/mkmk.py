# test of tests

from genice2.plugin import scan
import json

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
        if type(test) is dict:
            options = test["options"]
            test = test["args"]
        else:
            options = ""
        product = f"{ice}_{i}.gro"
        target = f"{ice}{test}"
        rules.append(f"{product}: ../../genice2/lattices/{ice}.py\n")
        rules.append(f"\t../../genice.x {target} {options} > $@\n")
        products.append(product)

print("all:",*products)
print("".join(rules))
