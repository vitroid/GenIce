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



def safe_import(category, name, **kwargs):
    logger = logging.getLogger()
    assert category in ("lattice", "format", "molecule", "loader")

    logger.warn((category, name, kwargs))


    usage = False
    if name[-1:] == "?":
        usage = True
        name = name[:-1]

    # kwargsが与えられた時はそちらを優先して利用する。
    # validate_argsがあればそれにそのまま渡す。
    # argparserがある場合には、後方互換性のために、kwargsを一旦文字列形式に戻し、それを渡す。
    # validate_argsとargparsertの両方を同時に指定するのは禁止。

    # name may contain arguments (obsolete)
    left = name.find("[")
    right = name.rfind("]")
    arg = ""
    if left < right:
        logger.warn("Deprecated way of specifying the module options.")
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

    module.logger = logger
    if usage:
        if "desc" in module.__dict__:
            logger.info("Usage for '{0}' plugin".format(name))
            print(module.desc["usage"])
            sys.exit(0)

    if category == "format":
        return module
    # 最終的には、この先のオプションチェックはすべて省く。
    logger.debug(category)
    # every plugin should have an argparser.
    if "argparser" in module.__dict__:
        logger.warn("Argparser is deprecated.")
        assert "validate_args" not in module.__dict__, "Use validate_args only."
        if arg == "" and len(kwargs) > 0:
            # backward compat; make arg from kwargs
            args = []
            for k, v in kwargs.items():
                if v is True:
                    args.append(k)
                else:
                    args.append("{0}={1}".format(k,v))
            arg = ":".join(args)
        module.arg    = arg
        # argparser for format plugin should be invoked in different place,
        # after the definition of a lattice object, as a hook.
        logger.info("Load {0} module '{1}', arguments [{2}]".format(category, name, arg))
        module.argparser(arg)
    elif "validate_args" in module.__dict__:
        logger.info("Load {0} module '{1}', arguments [{2}]".format(category, name, kwargs))
        invalid = module.validate_args(kwargs)
        if invalid:
            logger.warn("Invalid options: {0}".format(invalid))
    elif arg != "" or len(kwargs) == 0:
        logger.warn("Options are given but the module does not accept them.")
        if "usage" in module.__dict__:
            module.usage()
        else:
            if module.__doc__ is not None:
                for line in module.__doc__.splitlines():
                    logger.info("  "+line)
    return module
