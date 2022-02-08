# test of tests
import sys
from genice2.plugin import scan
import json

testmode = (len(sys.argv) > 1 and sys.argv[1] == "test")


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
        if testmode:
            rules.append(f"{product}.diff: {product} ../../genice2/lattices/{ice}.py\n")
            rules.append(f"\t../../genice.x {target} {options} | diff - $< \n") #> $@.diff\n")
            products.append(f"{product}.diff")
        else:
            rules.append(f"{product}: ../../genice2/lattices/{ice}.py\n")
            rules.append(f"\t../../genice.x {target} {options} > $@\n")
            products.append(product)

print("all:",*products)
print("".join(rules))
