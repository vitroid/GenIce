#!/usr/bin/env python
import sys
import os
# from genice2.tool import line_replacer
import distutils.core
from logging import getLogger, INFO, basicConfig
from jinja2 import Template
import json
from genice2.plugin import plugin_descriptors



def make_citations(r):
    if len(r) > 0:
        return " [" + ", ".join(sorted(r, key=lambda x: x.split()[-1])) + "]"
    return ""


def system_ices(markdown=True, citations=None):
    desc = plugin_descriptors("lattice", groups=["system"])
    documented, undocumented, refss = desc["system"]

    s = ""
    for description, ices in documented.items():
        citation = make_citations(refss[description])
        s += ", ".join(ices) + " | " + description + citation + "\n"
        if citations is not None:
            for ref in refss[description]:
                assert ref in citations, f"{ref} in {ices}"
    s += ", ".join(undocumented) + " | (Undocumented)\n"
    return s

def system_molecules(markdown=True, water=False, citations=None):
    desc = plugin_descriptors("molecule", water=water, groups=["system"])
    documented, undocumented, refss = desc["system"]

    s = ""
    for description, ices in documented.items():
        citation = make_citations(refss[description])
        s += ", ".join(ices) + " | " + description + citation + "\n"
        if citations is not None:
            for ref in refss[description]:
                assert ref in citations, f"{ref} in {ices}"
    s += ", ".join(undocumented) + " | (Undocumented)\n"
    return s



basicConfig(level=INFO, format="%(levelname)s %(message)s")
logger = getLogger()
logger.debug("Debug mode.")

with open("citations.json") as f:
    citations = json.load(f)

citationlist = [f"[{key}] {desc}" for key, doi, desc in citations]

def prefix(L, pre):
    return pre + ("\n"+pre).join(L) + "\n"

setup = distutils.core.run_setup("setup.py")

d = {
    "usage"   : prefix([x.rstrip() for x in os.popen("./genice.x -h").readlines()], "    "),
    "version" : setup.get_version(),
    "package" : setup.get_name(),
    "url"     : setup.get_url(),
    "genice"  : "[GenIce](https://github.com/vitroid/GenIce)",
    "requires": prefix(setup.install_requires, "* "),
    "ices"    : system_ices(), #citations=[key for key, doi, desc in citations]),
    "waters"  : system_molecules(water=True),
    "guests"  : system_molecules(water=False),
    "citationlist": prefix(citationlist, "* ")
}


t = Template(sys.stdin.read())
print(t.render(**d))
