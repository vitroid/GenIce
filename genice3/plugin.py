"""
Plugin handler.
"""

import glob
import importlib
import os
import re
import sys
import inspect
from collections import defaultdict
from functools import wraps
from logging import DEBUG, INFO, basicConfig, getLogger
from textwrap import fill
from typing import Any, Callable, Dict, Tuple

# import pkg_resources as pr

# from genice2.decorators import banner, timeit
from genice3.optionparser import parse_options

if sys.version_info < (3, 10):
    from importlib_metadata import entry_points
else:
    from importlib.metadata import entry_points


def parse_plugin_options(func: Callable) -> Callable:
    """
    プラグイン関数のオプション引数を自動的にパースするデコレータ

    このデコレータを使用すると、プラグイン関数に文字列のオプションを渡した場合、
    自動的に辞書形式にパースされます。

    使用方法:
        @parse_plugin_options
        def dump(genice: GenIce3, file: TextIOWrapper, options: str = ""):
            # optionsは自動的に辞書に変換されている
            guest = options.get("guest", {})
            ...

    または、キーワード引数として:
        @parse_plugin_options
        def dump(genice: GenIce3, file: TextIOWrapper, **options):
            # optionsは自動的に辞書に変換されている（options引数が文字列の場合）
            guest = options.get("guest", {})
            ...

    Args:
        func: デコレートするプラグイン関数

    Returns:
        ラップされた関数（オプション引数が自動的にパースされる）

    Examples:
        >>> @parse_plugin_options
        ... def dump(genice, file, options=""):
        ...     return options
        ...
        >>> dump(None, None, options="guest.A12=me,shift=(0.1,0.1,0.1)")
        {'guest': {'A12': 'me'}, 'shift': [0.1, 0.1, 0.1]}
    """
    sig = inspect.signature(func)

    @wraps(func)
    def wrapper(*args, **kwargs):
        # genice3.optionparserを遅延インポート（循環インポートを避けるため）
        from genice3.optionparser import parse_options

        # 関数の引数をバインド
        bound_args = sig.bind(*args, **kwargs)
        bound_args.apply_defaults()

        # 'options'という名前の引数が文字列の場合はパース
        if "options" in bound_args.arguments:
            options_value = bound_args.arguments["options"]
            if isinstance(options_value, str):
                # 文字列をパースして辞書に変換（空文字列の場合は空辞書）
                parsed_options = parse_options(options_value) if options_value else {}
                bound_args.arguments["options"] = parsed_options
                # kwargsにも反映
                if "options" in kwargs:
                    kwargs["options"] = parsed_options

        # 元の関数を呼び出し
        return func(*bound_args.args, **bound_args.kwargs)

    return wrapper


def scan(category):
    """
    Scan available plugins.
    """
    logger = getLogger()

    modules = {}
    desc = dict()
    iswater = dict()
    refs = dict()
    tests = dict()

    logger.info(f"\nPredefined {category}")
    module = importlib.import_module(f"genice3.{category}")
    mods = []
    for path in module.__path__:
        for mod in sorted(glob.glob(path + "/*.py")):
            mod = os.path.basename(mod)[:-3]
            if mod[:2] != "__":
                mods.append(mod)
    logger.info(mods)
    modules["system"] = mods

    for mod in modules["system"]:
        try:
            module = importlib.import_module(f"genice3.{category}.{mod}")
            if "desc" in module.__dict__:
                desc[mod] = module.desc["brief"]
                if "ref" in module.desc:
                    refs[mod] = module.desc["ref"]
                if "test" in module.desc:
                    tests[mod] = module.desc["test"]
            iswater[mod] = "water" in module.__dict__
        except BaseException:
            pass

    logger.info(f"Extra {category}")
    groupname = f"genice2_{category}"
    mods = []
    # for ep in pr.iter_entry_points(group=groupname):
    for ep in entry_points(group=groupname):
        mods.append(ep.name)
        try:
            module = ep.load()
            if "desc" in module.__dict__:
                desc[ep.name] = module.desc["brief"]
                if "ref" in module.desc:
                    refs[ep.name] = module.desc["ref"]
                if "test" in module.desc:
                    tests[mod] = module.desc["test"]
            iswater[ep.name] = "water" in module.__dict__
        except BaseException:
            pass
    logger.info(mods)
    modules["extra"] = mods

    logger.info(f"Local {category}")
    mods = [
        os.path.basename(mod)[:-3] for mod in sorted(glob.glob(f"./{category}/*.py"))
    ]
    logger.info(mods)
    for mod in mods:
        module = importlib.import_module(f"{category}.{mod}")
        if "desc" in module.__dict__:
            desc[mod] = module.desc["brief"]
            if "ref" in module.desc:
                refs[mod] = module.desc["ref"]
            if "test" in module.desc:
                tests[mod] = module.desc["test"]
        iswater[mod] = "water" in module.__dict__
    logger.info(mods)
    modules["local"] = mods
    modules["desc"] = desc
    modules["iswater"] = iswater
    modules["refs"] = refs
    modules["tests"] = tests

    return modules


def descriptions(category, width=72, water=False, groups=("system", "extra", "local")):
    """
    Show the list of available plugins in the category.

    Options:
      width=72      Width of the output.
      water=False   Pick up water molecules only (for molecule plugin).
    """
    titles = {
        "lattice": {
            "system": "1. Lattice structures served with GenIce",
            "extra": "2. Lattice structures served by external plugins",
            "local": "3. Lattice structures served locally",
            "title": "[Available lattice structures]",
        },
        "format": {
            "system": "1. Formatters served with GenIce",
            "extra": "2. Formatters served by external plugins",
            "local": "3. Formatters served locally",
            "title": "[Available formatters]",
        },
        "loader": {
            "system": "1. File types served with GenIce",
            "extra": "2. File types served by external eplugins",
            "local": "3. File types served locally",
            "title": "[Available input file types]",
        },
        "molecule": {
            "system": "1. Molecules served with GenIce",
            "extra": "2. Molecules served by external plugins",
            "local": "3. Molecules served locally",
            "title": "[Available molecules]",
        },
    }
    mods = scan(category)
    catalog = f" \n \n{titles[category]['title']}\n \n"
    desc = mods["desc"]
    iswater = mods["iswater"]
    for group in groups:
        desced = defaultdict(list)
        undesc = []
        for L in mods[group]:
            if category == "molecule":
                if L not in iswater:
                    iswater[L] = False
                if water and not iswater[L]:
                    continue
                if not water and iswater[L]:
                    continue
            if L in desc:
                desced[desc[L]].append(L)
            else:
                undesc.append(L)
        for dd in desced:
            desced[dd] = ", ".join(desced[dd])
        catalog += f"{titles[category][group]}\n \n"
        table = ""
        for dd in sorted(desced, key=lambda x: desced[x]):
            table += f"{desced[dd]}\t{dd}\n"
        if table == "":
            table = "(None)\n"
        table = [
            fill(
                line,
                width=width,
                drop_whitespace=False,
                expand_tabs=True,
                tabsize=16,
                subsequent_indent=" " * 16,
            )
            for line in table.splitlines()
        ]
        table = "\n".join(table) + "\n"
        undesc = " ".join(undesc)
        if undesc != "":
            undesc = "(Undocumented) " + undesc
        catalog += table + "----\n" + undesc + "\n \n \n"
    return catalog


def plugin_descriptors(category, water=False, groups=("system", "extra", "local")):
    """
    Show the list of available plugins in the category.

    Options:
      water=False   Pick up water molecules only (for molecule plugin).
    """
    mods = scan(category)
    catalog = dict()
    desc = mods["desc"]
    iswater = mods["iswater"]
    refs = mods["refs"]
    for group in groups:
        desced = defaultdict(list)
        undesc = []
        refss = defaultdict(set)
        for L in mods[group]:
            if category == "molecule":
                if L not in iswater:
                    iswater[L] = False
                if water and not iswater[L]:
                    continue
                if not water and iswater[L]:
                    continue
            if L in desc:
                # desc[L] is the brief description of the module
                # L is the name of module (name of ice)
                desced[desc[L]].append(L)
                if L in refs:
                    refss[desc[L]] |= set([label for key, label in refs[L].items()])
            else:
                undesc.append(L)
        catalog[group] = [desced, undesc, refss]
    return catalog


def audit_name(name):
    """
    Audit the mol name to avoid the access to external files
    """
    return re.match("^[A-Za-z0-9-_]+$", name) is not None


def import_extra(category, name):
    logger = getLogger()
    logger.info(f"Extra {category} plugin: {name}")
    groupname = f"genice2_{category}"
    module = None
    # for ep in pr.iter_entry_points(group=groupname):
    for ep in entry_points(group=groupname):
        logger.debug(f"    Entry point: {ep}")
        if ep.name == name:
            logger.debug(f"      Loading {name}...")
            module = ep.load()
    assert module is not None, f"Nonexistent or failed to load the module: {name}"
    return module


def safe_import(category, name):
    """
    Load a plugin.

    The plugins can exist either in the system, as a extra plugin, or in the
    local folder.

    category: The type of the plugin; "lattice", "format", "molecule", or "loader".
    name:     The name of the plugin.
    """
    logger = getLogger()
    assert category in ("exporter", "molecule", "unitcell")

    # single \? as a plugin name ==> show descriptions
    if name == "?":
        print(descriptions(category))
        sys.exit(0)

    usage = False
    if name[-1:] == "?":
        usage = True
        name = name[:-1]

    assert audit_name(name), f"Dubious {category} name: {name}"

    module = None
    fullname = f"{category}.{name}"
    logger.debug(f"Try to Load a local module: {fullname}")
    try:
        module = importlib.import_module(fullname)  # at ~/.genice
        logger.debug("Succeeded.")
    except ModuleNotFoundError:
        logger.debug(f"Module not found: {fullname}")
        module = None
    except ImportError as e:
        logger.error(f"Error importing module {fullname}: {str(e)}")
        raise
    if module is None:
        fullname = f"genice3.{category}.{name}"
        logger.debug(f"Try to load a system module: {fullname}")
        try:
            module = importlib.import_module(fullname)
            logger.debug("Succeeded.")
        except ModuleNotFoundError:
            logger.debug(f"Module not found: {fullname}")
            module = None
        except ImportError as e:
            logger.error(f"Error importing module {fullname}: {str(e)}")
            raise
    if module is None:
        logger.debug(f"Try to load an extra module: {fullname}")
        module = import_extra(category, name)
        logger.debug("Succeeded.")

    if usage:
        if "desc" in module.__dict__:
            logger.info(f"Usage for '{name}' plugin")
            print(module.desc["usage"])
            sys.exit(0)

    # exporterプラグインの場合、dump関数に自動的にデコレータを適用
    if category == "exporter" and hasattr(module, "dump"):
        module.dump = parse_plugin_options(module.dump)

    return module


def parse_plugin_specification(spec_str: str) -> tuple[str, Dict[str, Any]]:
    """
    commandline文字列をパースして、プラグイン名とパース済みオプション辞書を返す

    Args:
        commandline: プラグイン名、または "plugin[options]" 形式の文字列

    Returns:
        (plugin_name, parsed_options)
        - parsed_options はパース済みのオプション辞書

    Examples:
        >>> parse_exporter_spec("gromacs")
        ('gromacs', {})

        >>> parse_exporter_spec("gromacs[guest.A12=me,shift=(0.1,0.1,0.1)]")
        ('gromacs', {'guest': {'A12': 'me'}, 'shift': [0.1, 0.1, 0.1]})
    """
    # 括弧形式からオプション文字列を抽出
    if "[" in spec_str and "]" in spec_str:
        left = spec_str.find("[")
        right = spec_str.rfind("]")  # 最後の ] を取得
        plugin_name = spec_str[:left]
        option_string = spec_str[left + 1 : right]

        # オプション文字列をパース

        parsed_options = parse_options(option_string) if option_string else {}
        return plugin_name, parsed_options
    else:
        # プラグイン名のみの場合
        return spec_str, {}


def Lattice(name, **kwargs):
    """
    Shortcut for safe_import.
    """
    return safe_import("lattice", name).Lattice(**kwargs)


def Format(name, **kwargs):
    """
    Shortcut for safe_import.
    """
    return safe_import("format", name).Format(**kwargs)


def Molecule(name, **kwargs):
    """
    Shortcut for safe_import.
    """
    return safe_import("molecule", name).Molecule(**kwargs)


def Loader(name, **kwargs):
    """
    Shortcut for safe_import.
    """
    return safe_import("loader", name).Loader(**kwargs)


def Group(name, **kwargs):
    """
    Shortcut for safe_import.
    """
    return safe_import("group", name).Group(**kwargs)


if __name__ == "__main__":
    basicConfig(level=INFO)
    if len(sys.argv) == 1:
        cats = ("lattice", "format", "molecule", "loader")
    else:
        cats = sys.argv[1:]
    modules = {cat: scan(cat) for cat in cats}
    print(modules)
