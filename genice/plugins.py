#!/usr/bin/env python
from logging import getLogger, basicConfig, INFO, DEBUG
import pkg_resources as pr
import genice.lattices
import genice.molecules
import genice.formats
import genice.loaders
import glob
import os
import importlib
import sys

def scan(category):
    logger = getLogger()

    modules = {}
    desc = dict()
    modules["desc"] = desc

    logger.info("\nPredefined {0}s".format(category))
    module = importlib.import_module("genice.{0}s".format(category))
    mods = []
    for path in module.__path__:
        for mod in sorted(glob.glob(path+"/*.py")):
            mod = os.path.basename(mod)[:-3]
            if mod[:2] != "__":
                mods.append(mod)
    logger.info(mods)
    modules["system"] = mods

    for mod in modules["system"]:
        try:
            module = importlib.import_module("genice.{0}s.{1}".format(category, mod))
            if "desc" in module.__dict__:
                desc[mod] = module.desc["brief"]
        except:
            pass
    
    logger.info("Extra {0}s".format(category))
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
    logger.info(mods)
    modules["extra"] = mods

    logger.info("Local {0}s".format(category))
    mods = [os.path.basename(mod)[:-3] for mod in sorted(glob.glob("./{0}s/*.py".format(category)))]
    for mod in mods:
        module = importlib.import_module("{0}s.{1}".format(category, mod))
        if "desc" in module.__dict__:
            desc[mod] = module.desc["brief"]
    logger.info(mods)
    modules["local"] = mods

    return modules

if __name__ == "__main__":
    basicConfig(level=INFO)
    if len(sys.argv) == 1:
        cats = ("lattice", "format", "molecule", "loader")
    else:
        cats = argv[1:]
    modules = {cat: scan(cat) for cat in cats}
    print(modules)

    
