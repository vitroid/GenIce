#!/usr/bin/env python

# standard plugins
from logging import getLogger, basicConfig, INFO, DEBUG
import pkg_resources as pr
from textwrap import fill, wrap
import re
from collections import defaultdict
import glob
import os
import importlib
import sys

# public libs
import genice.lattices
import genice.molecules
import genice.formats
import genice.loaders

def scan(category):
    logger = getLogger()

    modules = {}
    desc = dict()
    modules["desc"] = desc

    logger.debug("Predefined {0}s".format(category))
    module = importlib.import_module("genice.{0}s".format(category))
    mods = []
    for path in module.__path__:
        for mod in sorted(glob.glob(path+"/*.py")):
            mod = os.path.basename(mod)[:-3]
            if mod[:2] != "__":
                mods.append(mod)
    logger.debug(mods)
    modules["system"] = mods

    for mod in modules["system"]:
        try:
            module = importlib.import_module("genice.{0}s.{1}".format(category, mod))
            if "desc" in module.__dict__:
                desc[mod] = module.desc["brief"]
        except:
            pass
    
    logger.debug("Extra {0}s".format(category))
    groupname = 'genice_{0}'.format(category)
    mods = []
    for ep in pr.iter_entry_points(group=groupname):
        label, m = str(ep).split("=")
        mods.append(label)
        try:
            module = ep.load()
            if "desc" in module.__dict__:
                desc[label] = module.desc["brief"]
        except:
            pass
    logger.debug(mods)
    modules["extra"] = mods

    logger.debug("Local {0}s".format(category))
    mods = [os.path.basename(mod)[:-3] for mod in sorted(glob.glob("./{0}s/*.py".format(category)))]
    for mod in mods:
        module = importlib.import_module("{0}s.{1}".format(category, mod))
        if "desc" in module.__dict__:
            desc[mod] = module.desc["brief"]
    logger.debug(mods)
    modules["local"] = mods

    return modules




def descriptions(category, width=72):
    titles={ "lattice": {"system": "1. Lattice structures served with GenIce",
                         "extra":  "2. Lattice structures served by external plugins",
                         "local":  "3. Lattice structures served locally",
                         "title":  "[Available lattice structures]"},
             "format": {"system": "1. Formatters served with GenIce",
                        "extra":  "2. Formatters served by external plugins",
                        "local":  "3. Formatters served locally",
                        "title":  "[Available formatters]"},
             "loader": {"system": "1. File types served with GenIce",
                        "extra":  "2. File types served by external eplugins",
                        "local":  "3. File types served locally",
                        "title":  "[Available input file types]"},
             "molecule": {"system": "1. Molecules served with GenIce",
                        "extra":  "2. Molecules served by extternal plugins",
                        "local":  "3. Molecules served locally",
                        "title":  "[Available molecules]"},
             }
    mods = scan(category)
    catalog = " \n \n{0}\n \n".format(titles[category]["title"])
    desc = mods["desc"]
    for group in ("system", "extra", "local"):
        desced = defaultdict(list)
        undesc = []
        for L in mods[group]:
            if L in desc:
                desced[desc[L]].append(L)
            else:
                undesc.append(L)
        for dd in desced:
            desced[dd] = ", ".join(desced[dd])
        catalog += "{0}\n \n".format(titles[category][group])
        table = ""
        for dd in sorted(desced, key=lambda x: desced[x]):
            table += "{0}\t{1}\n".format(desced[dd], dd)
        if table == "":
            table = "(None)\n"
        table = [fill(line, width=width, drop_whitespace=False, expand_tabs=True, tabsize=16, subsequent_indent=" "*16) for line in table.splitlines()]
        table = "\n".join(table)+"\n"
        undesc = " ".join(undesc)
        if undesc != "":
            undesc = "(Undocumented) " + undesc
        catalog += table + "----\n" + undesc + "\n \n \n"
    return catalog


def audit_name(name):
    """
    Audit the mol name to avoid the access to external files
    """
    return re.match('^[A-Za-z0-9-_]+$', name) is not None


def import_extra(category, name):
    logger = getLogger()
    logger.info("Extra {0} plugin: {1}".format(category,name))
    groupname = 'genice_{0}'.format(category)
    for ep in pr.iter_entry_points(group=groupname):
        logger.debug("    Entry point: {0}".format(ep))
        if ep.name == name:
            logger.debug("      Loading {0}...".format(name))
            module = ep.load()
    return module



def safe_import(category, name):
    logger = getLogger()
    assert category in ("lattice", "format", "molecule", "loader")

    # single \? as a plugin name ==> show descriptions
    if name == "?":
        print(descriptions(category))
        sys.exit(0)

    usage = False
    if name[-1:] == "?":
        usage = True
        name = name[:-1]

        
    # name may contain arguments
    left = name.find("[")
    right = name.rfind("]")
    arg = ""
    if left < right:
        arg = name[left + 1:right]
        name = name[:left]
    assert audit_name(name), "Dubious {0} name: {1}".format(category, name)

#    if category == 'format':
#        logger.info("Load {0} module '{1}', arguments [{2}]".format(category, name, arg))
#        return import_format_plugin(category, name), arg

    module = None
    try:
        module = importlib.import_module(category + "s." + name)  # at ~/.genice
    except ImportError as e:
        pass
    if module is None:
        fullname = "genice." + category + "s." + name
        logger.debug("Load module: {0}".format(fullname))
        try:
            module = importlib.import_module(fullname)
        except:
            pass
    if module is None:
        module = import_extra(category, name)
        
    logger.info("Load {0} module '{1}', arguments [{2}]".format(category, name, arg))
    module.arg    = arg
    if usage:
        if "desc" in module.__dict__:
            logger.info("Usage for '{0}' plugin".format(name))
            print(module.desc["usage"])
            sys.exit(0)

    if category == "format":
        return module
    logger.debug(category)
    # every plugin should have an argparser.
    if "argparser" in module.__dict__:
        # argparser for format plugin should be invoked in different place,
        # after the definition of a lattice object, as a hook.
        module.argparser(arg)
    elif arg != "":
        logger.info("Arguments are given but the module does not accept them.")
        if "usage" in module.__dict__:
            module.usage()
        else:
            for line in module.__doc__.splitlines():
                logger.info("  "+line)
    return module




if __name__ == "__main__":
    basicConfig(level=INFO)
    if len(sys.argv) == 1:
        cats = ("lattice", "format", "molecule", "loader")
    else:
        cats = argv[1:]
    modules = {cat: scan(cat) for cat in cats}
    print(modules)

    
