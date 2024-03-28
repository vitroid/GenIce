from genice2.plugin import scan, safe_import


def test_lattice_all():
    def test_required_attributes(ice):
        # assert hasattr(ice, "density")
        assert hasattr(ice, "waters")
        assert hasattr(ice, "coord")
        assert hasattr(ice, "cell")

    def test_available_attributes(ice):
        for member in dir(ice):
            if not member.startswith("__"):
                assert member in (
                    "waters",
                    "coord",
                    "cell",
                    "fixed",
                    "density",
                    "bondlen",
                    "pairs",
                    "cagepos",
                    "cagetype",
                    "cages",
                )

    category = "lattice"
    plugins = scan(category)

    for plugin_name in plugins["system"]:
        plugin = safe_import(category, plugin_name)
        # if plugin_name in plugins["desc"]:
        #     # it has special test suite
        #     print(plugin_name, plugins["desc"][plugin_name])

        # pluginのオプション
        test_kwargs = []

        # if the plugin has its own test cases:
        if plugin_name in plugins["tests"]:
            for testcase in plugins["tests"][plugin_name]:
                if "args" in testcase:
                    # print(plugin_name, plugins["tests"][plugin_name], " test")
                    if type(testcase["args"]) is str:
                        if testcase["args"] != "":
                            test_kwargs.append({testcase["args"]: True})
                    else:
                        assert type(testcase["args"]) is dict
                        test_kwargs.append(testcase["args"])
        if len(test_kwargs) == 0:
            test_kwargs = [{}]
        else:
            print(plugin_name, test_kwargs)
        for kwarg in test_kwargs:
            ice = plugin.Lattice(**kwarg)
            test_required_attributes(ice)
            test_available_attributes(ice)


if __name__ == "__main__":
    test_lattice_all()
