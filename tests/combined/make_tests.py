from genice2.plugin import scan, safe_import
from logging import getLogger, basicConfig, DEBUG
import sys


def formats():
    category = "format"
    plugins = scan(category)

    for plugin_name in plugins["system"]:
        plugin = safe_import(category, plugin_name)
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

    category = "lattice"
    plugins = scan(category)

    for plugin_name in plugins["system"]:
        plugin = safe_import(category, plugin_name)
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
                    yield plugin_name, testcase
                    yieldcount += 1
                elif "options" in testcase:
                    testcase["args"] = {}
                    yield plugin_name, testcase
                    yieldcount += 1
        if yieldcount == 0:
            yield plugin_name, {"args": {}}


def waters():

    category = "molecule"
    plugins = scan(category)

    for plugin_name in plugins["system"]:
        plugin = safe_import(category, plugin_name)
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
                    if "water" in dir(plugin) and plugin.water == 1:
                        yield plugin_name, testcase["args"]
                        yieldcount += 1
        if yieldcount == 0:
            if "water" in dir(plugin) and plugin.water == 1:
                yield plugin_name, {}


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
    for i, (water, lattice, format) in enumerate(
        zip(startover(waters), lattices(), startover(formats))
    ):
        args = ["$(GENICE)"]
        args += [lattice[0] + make_plugin_args(lattice[1]["args"])]
        args += ["-w", water[0] + make_plugin_args(water[1])]
        args += ["-f ", format[0] + make_plugin_args(format[1])]
        if "options" in lattice[1]:
            args.append(lattice[1]["options"])
        if lattice[0] == "8":
            print(lattice)

        yield f"_test{i}", " ".join(args)


if __name__ == "__main__":
    # basicConfig(level=DEBUG)
    prefix = sys.argv[1]

    all = []

    s = ""
    for filename, rule in makefile_rules():
        s += prefix + "/" + filename + ":\n\t" + rule + " > $@\n"
        all.append(prefix + "/" + filename)

    print("GENICE=../../genice.x")
    print("all: " + " ".join(all))
    print(s)
