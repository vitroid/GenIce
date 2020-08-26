from logging import getLogger

def line_replacer(line, d):
    logger = getLogger()
    s = ""
    for tag in d:
        loc = line.find(tag)
        if loc >= 0:
            logger.debug("From {0} by {1}.".format(tag, d[tag]))
            replacement = d[tag].splitlines()
            if len(replacement) == 1:
                s = line.replace(tag, replacement[0])
            else:
                indent = line[:loc]
                for newline in replacement:
                    s += indent + newline + "\n"
            return s
    return line



def plugin_option_parser(s):
    """
    Separate the plugin name and options
    """

    left = s.find("[")
    right = s.find("]")
    if 0 < left < len(s) and 0 < right < len(s) and left < right:
        args = s[left+1:right]
        name = s[:left]
    else:
        return s, {}

    kwargs = dict()
    for elem in args.split(":"):
        if "=" in elem:
            k, v = elem.split("=", 2)
            kwargs[k] = v
        else:
            kwargs[elem] = True
    return name, kwargs
