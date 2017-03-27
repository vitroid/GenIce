import re
import importlib
import logging

def audit_name(name):
    """
    Audit the mol name to avoid the access to external files
    """
    return re.match('^[A-Za-z0-9-_]+$', name) is not None



def safe_import(category, name):
    logger = logging.getLogger()
    assert category in ("lattice", "format", "molecule")
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
    return module
