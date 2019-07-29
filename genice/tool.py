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

