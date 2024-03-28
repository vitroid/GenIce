from genice2.plugin import scan, safe_import
from logging import getLogger, basicConfig, DEBUG


def test_molecule_all():
    def test_required_attributes(molecule):
        # assert hasattr(ice, "density")
        assert hasattr(molecule, "sites_")
        assert hasattr(molecule, "labels_")
        assert hasattr(molecule, "name_")

    def test_available_attributes(molecule):
        for member in dir(molecule):
            if not member.startswith("__"):
                assert member in (
                    "sites_",
                    "labels_",
                    "name_",
                    "atoms_",
                    "get",
                ), f"The plaugin {plugin_name} defines an extra variable {member}."

    def test_water(molecule):
        for label, atom in zip(molecule.labels_[:3], "OHH"):
            assert label[0] == atom, f"First three atoms of water must be O, H, and H."

    category = "molecule"
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
            molecule = plugin.Molecule(**kwarg)
            test_required_attributes(molecule)
            test_available_attributes(molecule)
        if "water" in dir(plugin) and plugin.water == 1:
            test_water(molecule)


if __name__ == "__main__":
    basicConfig(level=DEBUG)
    test_molecule_all()
