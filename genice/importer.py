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



def import_format_extra(name):
    logger = logging.getLogger()
    logger.info("Extra plugin: {0}".format(name))
    hooks = dict()
    for i in (0,1,2,3,4,5,6,7):
        logger.debug("  Looking for hook{0}...".format(i))
        groupname = 'genice_format_hook{0}'.format(i)
        for ep in pr.iter_entry_points(group=groupname):
            logger.debug("    Entry point for hook{1}: {0}".format(ep, i))
            if ep.name == name:
                logger.debug("      Loading {0}...".format(name))
                hooks[i] = ep.load()
                logger.info("      Loaded {1} Hook{0}.".format(i,name))
    if len(hooks) == 0:
        logger.error("No extra plugin named '{0}'.".format(name))
        sys.exit(1)
    return hooks


def import_format_plugin(category, name):
    logger = logging.getLogger()
    module = None
    try:
        module     = importlib.import_module(category+"s."+name) #at ~/.genice
    except ImportError as e:
        pass
    if module is None:
        fullname = "genice."+category+"s."+name
        logger.debug("Load module: {0}".format(fullname))
        try:
            module     = importlib.import_module(fullname)
        except ImportError as e:
            pass
    if module is None:
        # load extras
        hooks = import_format_extra(name)
    else:
        # load user-defined.
        hooks = module.hooks
    return hooks
    


def safe_import(category, name):
    logger = logging.getLogger()
    assert category in ("lattice", "format", "molecule", "loader")
    # name may contain arguments
    left = name.find("[")
    right = name.rfind("]")
    arg = ""
    if left < right:
        arg = name[left+1:right]
        name = name[:left]
    assert audit_name(name), "Dubious {0} name: {1}".format(category, name)

    if category == 'format':
        logger.info("Load {0} module '{1}', arguments [{2}]".format(category, name, arg))
        return import_format_plugin(category, name), arg

    module = None
    try:
        module     = importlib.import_module(category+"s."+name) #at ~/.genice
    except ImportError as e:
        pass
    if module is None:
        fullname = "genice."+category+"s."+name
        logger.debug("Load module: {0}".format(fullname))
        module     = importlib.import_module(fullname)
    logger.info("Load {0} module '{1}', arguments [{2}]".format(category, name, arg))
    if arg != "":
        if "argparser" in module.__dict__:
            module.argparser(arg)
        else:
            logger.info("Arguments are given but the module does not accept them.")
    elif "usage" in module.__dict__:
        module.usage()
    if category == "lattice":
        module.lattice_type = name
    return module
