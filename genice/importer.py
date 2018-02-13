import re
import importlib
import logging

#
# Accept options parenthesized after the plugin name.
#


def audit_name(name):
    """
    Audit the mol name to avoid the access to external files
    """
    return re.match('^[A-Za-z0-9-_]+$', name) is not None



def safe_import(category, name):
    logger = logging.getLogger()
    assert category in ("lattice", "format", "molecule")
    # name may contain arguments
    left = name.find("[")
    right = name.rfind("]")
    arg = ""
    if left < right:
        arg = name[left+1:right]
        name = name[:left]
    assert audit_name(name), "Dubious {0} name: {1}".format(category, name)
    module = None
    try:
        module     = importlib.import_module(category+"s."+name) #at ~/.genice
    except ImportError as e:
        pass
    if module is None:
        fullname = "genice."+category+"s."+name
        logger.debug("Load module: {0}".format(fullname))
        module     = importlib.import_module(fullname)
    logger.info("Load {0} module {1} arguments [{2}]".format(category, name, arg))
    if arg != "":
        if "argparser" in module.__dict__:
            module.argparser(arg)
        else:
            logger.info("Arguments are given but the module does not accept them.")
    elif "usage" in module.__dict__:
        module.usage()
    return module
