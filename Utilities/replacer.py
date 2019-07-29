#!/usr/bin/env python
import sys
import os
from genice.tool import line_replacer
import distutils.core

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
