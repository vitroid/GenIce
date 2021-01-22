#!/usr/bin/env python
import sys
import os
# from genice2.tool import line_replacer
import distutils.core
from logging import getLogger
from jinja2 import Template
import json
from genice2.plugin import plugin_descriptors



def system_ices(markdown=True):
    desc = plugin_descriptors("lattice", groups=["system"])
    documented, undocumented = desc["system"]

    s = ""
    for description, ices in documented.items():
        s += ", ".join(ices) + " | " + description + "\n"
    s += ", ".join(undocumented) + " | (Undocumented)\n"
    return s




with open("citations.json") as f:
    citations = json.load(f)

citationlist = [f"[{key}] {desc}" for key, doi, desc in citations]

def prefix(L, pre):
    return pre + ("\n"+pre).join(L) + "\n"

setup = distutils.core.run_setup("setup.py")

d = {
    "usage"   : prefix(os.popen("./genice.x -h").readlines(), "    "),
    "version" : setup.get_version(),
    "package" : setup.get_name(),
    "url"     : setup.get_url(),
    "genice"  : "[GenIce](https://github.com/vitroid/GenIce)",
    "requires": prefix(setup.install_requires, "* "),
    "ices"    : system_ices(),
    "citationlist": prefix(citationlist, "* ")
}


t = Template(sys.stdin.read())
print(t.render(**d))
