#!/usr/bin/env python
import sys
import os
# from genice2.tool import line_replacer
import distutils.core
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

setup = distutils.core.run_setup("setup.py")

d = {
    "%%usage%%"   : "".join(os.popen("./genice.x -h").readlines()),
    "%%usage_analice%%" : "".join(os.popen("./analice.x -h").readlines()),
    "%%version%%" : setup.get_version(),
    "%%package%%" : setup.get_name(),
    "%%url%%"     : setup.get_url(),
    "%%genice%%"  : "[GenIce](https://github.com/vitroid/GenIce)",
    "%%requires%%": "\n".join(setup.install_requires),
}


for line in sys.stdin:
    print(line_replacer(line, d), end="")
