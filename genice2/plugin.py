#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Plugin handler.
"""

import glob
import importlib
import os
import re
import sys
from collections import defaultdict
from logging import DEBUG, INFO, basicConfig, getLogger
from textwrap import fill

import pkg_resources as pr

from genice2.decorators import banner, timeit


@timeit
@banner
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

    logger.info(f"\nPredefined {category}s")
    module = importlib.import_module(f"genice2.{category}s")
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
            module = importlib.import_module(
               f"genice2.{category}s.{mod}")
            if "desc" in module.__dict__:
                desc[mod] = module.desc["brief"]
                if "ref" in module.desc:
                    refs[mod] = module.desc["ref"]
                if "test" in module.desc:
                    tests[mod] = module.desc["test"]
            iswater[mod] = "water" in module.__dict__
        except BaseException:
            pass

    logger.info(f"Extra {category}s")
    groupname = f'genice2_{category}'
    mods = []
    for ep in pr.iter_entry_points(group=groupname):
        label, m = str(ep).split("=")
        mods.append(label)
        try:
            module = ep.load()
            if "desc" in module.__dict__:
                desc[label] = module.desc["brief"]
                if "ref" in module.desc:
                    refs[label] = module.desc["ref"]
                if "test" in module.desc:
                    tests[mod] = module.desc["test"]
            iswater[label] = "water" in module.__dict__
        except BaseException:
            pass
    logger.info(mods)
    modules["extra"] = mods

    logger.info(f"Local {category}s")
    mods = [os.path.basename(mod)[:-3]
            for mod in sorted(glob.glob(f"./{category}s/*.py"))]
    logger.info(mods)
    for mod in mods:
        module = importlib.import_module(f"{category}s.{mod}")
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


def descriptions(
    category,
    width=72,
    water=False,
    groups=(
        "system",
        "extra",
        "local")):
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
            "title": "[Available lattice structures]"},
        "format": {
            "system": "1. Formatters served with GenIce",
            "extra": "2. Formatters served by external plugins",
            "local": "3. Formatters served locally",
            "title": "[Available formatters]"},
        "loader": {
            "system": "1. File types served with GenIce",
            "extra": "2. File types served by external eplugins",
            "local": "3. File types served locally",
            "title": "[Available input file types]"},
        "molecule": {
            "system": "1. Molecules served with GenIce",
            "extra": "2. Molecules served by external plugins",
            "local": "3. Molecules served locally",
            "title": "[Available molecules]"},
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
                subsequent_indent=" " *
                16) for line in table.splitlines()]
        table = "\n".join(table) + "\n"
        undesc = " ".join(undesc)
        if undesc != "":
            undesc = "(Undocumented) " + undesc
        catalog += table + "----\n" + undesc + "\n \n \n"
    return catalog


def plugin_descriptors(
    category,
    water=False,
    groups=(
        "system",
        "extra",
        "local")):
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
                    refss[desc[L]] |= set(
                        [label for key, label in refs[L].items()])
            else:
                undesc.append(L)
        catalog[group] = [desced, undesc, refss]
    return catalog


def audit_name(name):
    """
    Audit the mol name to avoid the access to external files
    """
    return re.match('^[A-Za-z0-9-_]+$', name) is not None


def import_extra(category, name):
    logger = getLogger()
    logger.info(f"Extra {category} plugin: {name}")
    groupname = f'genice2_{category}'
    module = None
    for ep in pr.iter_entry_points(group=groupname):
        logger.debug(f"    Entry point: {ep}")
        if ep.name == name:
            logger.debug(f"      Loading {name}...")
            module = ep.load()
    if module is None:
        logger.error(f"Nonexistent module: {name}")
        sys.exit(1)
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
    assert category in ("lattice", "format", "molecule", "loader", "group")

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
    try:
        module = importlib.import_module(
            category + "s." + name)  # at ~/.genice
    except ImportError as e:
        pass
    if module is None:
        fullname = "genice2." + category + "s." + name
        logger.debug(f"Load module: {fullname}")
        try:
            module = importlib.import_module(fullname)
        except BaseException:
            pass
    if module is None:
        module = import_extra(category, name)

    if usage:
        if "desc" in module.__dict__:
            logger.info(f"Usage for '{name}' plugin")
            print(module.desc["usage"])
            sys.exit(0)

    return module


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
