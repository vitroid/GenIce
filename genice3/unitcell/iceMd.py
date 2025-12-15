#!/usr/bin/python
# coding: utf-8

from genice3.unitcell.iceM import UnitCell as UnitCellM
import networkx as nx

desc = {
    "ref": {"Md": "Mochizuki 2024"},
    "usage": "No options available.",
    "brief": "A hydrogen-disordered counterpart of ice M.",
    "test": ({"options": "--depol=none"},),
}


class UnitCell(UnitCellM):
    def __init__(self):
        super().__init__()
        self.fixed = nx.DiGraph()
