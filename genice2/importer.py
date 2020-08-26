import re
import importlib
import logging
import pkg_resources as pr
import sys
#
# Accept options parenthesized after the plugin name.
#


def audit_name(name):
    """
    Audit the mol name to avoid the access to external files
    """
    return re.match('^[A-Za-z0-9-_]+$', name) is not None


def import_extra(category, name):
    logger = logging.getLogger()
    logger.info("Extra {0} plugin: {1}".format(category,name))
    groupname = 'genice_{0}'.format(category)
    module = None
    for ep in pr.iter_entry_points(group=groupname):
        logger.debug("    Entry point: {0}".format(ep))
        if ep.name == name:
            logger.debug("      Loading {0}...".format(name))
            module = ep.load()
    if module is None:
        logger.error("Nonexistent module: {0}".format(name))
        sys.exit(1)
    return module



def safe_import(category, name):
    logger = logging.getLogger()
    assert category in ("lattice", "format", "molecule", "loader")

    usage = False
    if name[-1:] == "?":
        usage = True
        name = name[:-1]

    assert audit_name(name), "Dubious {0} name: {1}".format(category, name)

    module = None
    try:
        module = importlib.import_module(category + "s." + name)  # at ~/.genice
    except ImportError as e:
        pass
    if module is None:
        fullname = "genice2." + category + "s." + name
        logger.debug("Load module: {0}".format(fullname))
        try:
            module = importlib.import_module(fullname)
        except:
            pass
    if module is None:
        module = import_extra(category, name)

    module.logger = logger
    if usage:
        if "desc" in module.__dict__:
            logger.info("Usage for '{0}' plugin".format(name))
            print(module.desc["usage"])
            sys.exit(0)

    return module
