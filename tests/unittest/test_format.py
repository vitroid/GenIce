from genice2.plugin import scan, safe_import
from logging import getLogger, basicConfig, DEBUG


def test_format_all():
    def test_required_attributes(format):
        # assert hasattr(ice, "density")
        assert hasattr(format, "hooks")

    def test_hook(format):
        for k, v in format.hooks().items():
            assert k in (1, 2, 3, 4, 5, 6, 7)
            assert callable(v)

    category = "format"
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
            format = plugin.Format(**kwarg)
            test_required_attributes(format)
            test_hook(format)


if __name__ == "__main__":
    basicConfig(level=DEBUG)
    test_format_all()
