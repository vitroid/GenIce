from genice2.plugin import scan, safe_import
from logging import getLogger, basicConfig, DEBUG
import sys
import random
import itertools as it


def formats():
    category = "format"
    plugins = scan(category)

    random.shuffle(plugins["system"])

    for plugin_name in plugins["system"]:
        # plugin = safe_import(category, plugin_name)
        # if plugin_name in plugins["desc"]:
        #     # it has special test suite
        #     print(plugin_name, plugins["desc"][plugin_name])

        yieldcount = 0

        # if the plugin has its own test cases:
        if plugin_name in plugins["tests"]:
            for testcase in plugins["tests"][plugin_name]:
                if "args" in testcase:
                    content = testcase["args"]
                    if type(content) is str:
                        if content != "":
                            testcase["args"] = {content: True}
                    assert type(testcase["args"]) is dict, plugin_name
                    yield plugin_name, testcase["args"]
                    yieldcount += 1
        if yieldcount == 0:
            yield plugin_name, {}


def lattices():

    def is_fixed(plugin, lattice_args={}):
        """if the lattice module has member "fixed" """
        lattice = plugin.Lattice(**lattice_args)
        # if the lattice module has member "fixed",
        return "fixed" in dir(lattice)

    def has_cages(plugin, lattice_args={}):
        """if the lattice module has member 'fixed'"""
        lattice = plugin.Lattice(**lattice_args)
        # if plugin_name == "engel29":
        #     assert False, dir(lattice)
        return "cages" in dir(lattice)

    category = "lattice"
    plugins = scan(category)

    random.shuffle(plugins["system"])

    for plugin_name in plugins["system"]:
        if plugin_name in ("iceR",):
            continue
        lattice = safe_import(category, plugin_name)
        # if plugin_name in plugins["desc"]:
        #     # it has special test suite
        #     print(plugin_name, plugins["desc"][plugin_name])

        # if the plugin has its own test cases:
        yieldcount = 0
        if plugin_name in plugins["tests"]:
            for testcase in plugins["tests"][plugin_name]:
                if "args" in testcase:
                    content = testcase["args"]
                    if type(content) is str:
                        if content != "":
                            testcase["args"] = {content: True}
                    assert type(testcase["args"]) is dict, plugin_name
                    if is_fixed(lattice, testcase["args"]):
                        # it is informed
                        testcase["fixed"] = True
                    if has_cages(lattice, testcase["args"]):
                        # it is informed
                        testcase["cages"] = True
                    yield plugin_name, testcase
                    yieldcount += 1
                elif "options" in testcase:
                    testcase["args"] = {}
                    if is_fixed(lattice):
                        # it is informed
                        testcase["fixed"] = True
                    if has_cages(lattice):
                        # it is informed
                        testcase["cages"] = True
                    yield plugin_name, testcase
                    yieldcount += 1
        if yieldcount == 0:
            testcase = {f"args": {}}
            if is_fixed(lattice):
                # it is informed
                testcase["fixed"] = True
            if has_cages(lattice):
                # it is informed
                testcase["cages"] = True
            yield plugin_name, testcase


def waters():

    category = "molecule"
    plugins = scan(category)

    random.shuffle(plugins["system"])

    for plugin_name in plugins["system"]:
        molecule = safe_import(category, plugin_name)
        # if plugin_name in plugins["desc"]:
        #     # it has special test suite
        #     print(plugin_name, plugins["desc"][plugin_name])

        yieldcount = 0

        # if the plugin has its own test cases:
        if plugin_name in plugins["tests"]:
            for testcase in plugins["tests"][plugin_name]:
                if "args" in testcase:
                    content = testcase["args"]
                    if type(content) is str:
                        if content != "":
                            testcase["args"] = {content: True}
                    assert type(testcase["args"]) is dict, plugin_name
                    if "water" in dir(molecule) and molecule.water == 1:
                        yield plugin_name, testcase["args"]
                        yieldcount += 1
        if yieldcount == 0:
            if "water" in dir(molecule) and molecule.water == 1:
                yield plugin_name, {}


general_options = [
    ["", "--rep 2 2 2", "--rep 3 3 3", "--reshape 1,1,1,1,-1,0,1,1,-2"],
    #  --assess_cages
    #  --Guest
    #  --guest
    ["", "--add_noise 2"],
    ["", "--dens 1.0"],
    ["", "--shift 2 2 2"],
]

# options for hydrogen-disordered ices only
hdi_options = ["", "--cation 0=NH4 --anion 5=F"]

hyd_options = ["--guest 12=g12*0.5+g14*0.3"]


def options():
    ops = [" ".join(x) for x in it.product(*general_options)]
    random.shuffle(ops)
    for op in ops:
        yield op


def startover(iter):
    while True:
        yield from iter()


def make_plugin_args(args):
    if len(args) == 0:
        return ""
    s = []
    for k, v in args.items():
        if v == True:
            s.append(f"{k}")
        else:
            s.append(f"{k}={v}")
    return "[" + ":".join(s) + "]"


def makefile_rules():
    for i, (water, (lattice_name, lattice_args), format, ops) in enumerate(
        zip(startover(waters), lattices(), startover(formats), startover(options))
    ):
        args = ["$(GENICE)"]
        args += [lattice_name + make_plugin_args(lattice_args["args"])]
        args += ["-w", water[0] + make_plugin_args(water[1])]
        args += ["-f ", format[0] + make_plugin_args(format[1])]
        args += [ops]
        if "fixed" not in lattice_args:
            args.append(random.choice(hdi_options))
        if "cages" in lattice_args:
            args.append(random.choice(hyd_options))
        if "options" in lattice_args:
            args.append(lattice_args["options"])

        yield f"_test{i}", " ".join(args)


if __name__ == "__main__":
    # basicConfig(level=DEBUG)
    seed = int(sys.argv[1])
    random.seed(seed)

    all = []

    s = ""
    for filename, rule in makefile_rules():
        s += filename + ":\n\t" + rule + " > $@\n"
        all.append(filename)

    print("GENICE=../../../genice.x")
    print("all: " + " ".join(all))
    print(s)
