from genice2.plugin import scan, safe_import
from logging import getLogger, basicConfig, DEBUG


input_files = {
    "ar3r": "3.ar3r",
    "exyz": "5.exyz",
    "gro": "7.gro",
    "mdv": "CS1.mdv",
    "mdva": "CS2.mdva",
    "nx3a": "HS1.nx3a",
}


def test_loader_all():

    category = "loader"
    plugins = scan(category)

    for plugin_name in plugins["system"]:
        plugin = safe_import(category, plugin_name)
        # if plugin_name in plugins["desc"]:
        #     # it has special test suite
        #     print(plugin_name, plugins["desc"][plugin_name])

        if plugin_name == "ar3a":
            kwargs = {}
        else:
            kwargs = {"oname": "OW", "hname": "HW*"}
        with open(input_files[plugin_name]) as file:
            molecules = plugin.load_iter(file, **kwargs)


if __name__ == "__main__":
    basicConfig(level=DEBUG)
    test_loader_all()
