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
    for ep in pr.iter_entry_points(group=groupname):
        logger.debug("    Entry point: {0}".format(ep))
        if ep.name == name:
            logger.debug("      Loading {0}...".format(name))
            module = ep.load()
    return module



def safe_import(category, name):
    logger = logging.getLogger()
    assert category in ("lattice", "format", "molecule", "loader")
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
    module.logger = logger
    module.arg    = arg

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
