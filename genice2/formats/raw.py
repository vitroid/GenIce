# coding: utf-8

desc={"ref": {},
      "brief": "Raw data.",
      "usage": """
It is convenient to select the raw data format when you call GenIce from
Jupyter and Google Colab environment. Raw data is not available in the command line.
"""}


from collections import defaultdict
from logging import getLogger
import numpy as np
from genice2.decorators import timeit, banner
import genice2.formats


class Format(genice2.formats.Format):
    def __init__(self, **kwargs):
        unknown = dict()
        self.stages = []
        for k, v in kwargs.items():
            if k == "stage":
                self.stages = v
            else:
                unknown[k] = v
        super().__init__(**unknown)
        self.output = dict() # a bucket of data


    def hooks(self):
        hook_funcs = {
            1:self.Hook1,
            2:self.Hook2,
            3:self.Hook3,
            4:self.Hook4,
            5:self.Hook5,
            6:self.Hook6,
            7:self.Hook7,
            }
        active = dict()
        for i in self.stages:
            if i in hook_funcs:
                active[i] = hook_funcs[i]
        return active


    @timeit
    @banner
    def Hook1(self, ice):
        "Replicate water molecules to make a repeated cell."
        logger = getLogger()
        self.output["reppositions"] = ice.reppositions
        self.output["repcell"]      = ice.repcell.mat
        self.output["repcagetype"]  = ice.repcagetype
        self.output["repcagepos"]   = ice.repcagepos
        self.output["cagetypes"]    = ice.cagetypes


    @timeit
    @banner
    def Hook2(self, ice):
        "Make a random graph and replicate."
        logger = getLogger()
        self.output["dopants"]      = ice.dopants
        self.output["groups"]       = ice.groups
        self.output["filled_cages"] = ice.filled_cages
        self.output["graph"]        = ice.graph


    @timeit
    @banner
    def Hook3(self, ice):
        "Make a true ice graph."
        logger = getLogger()
        self.output["graph"]        = ice.graph


    @timeit
    @banner
    def Hook4(self, ice):
        "Depolarize."
        logger = getLogger()
        self.output["spacegraph"]   = ice.spacegraph


    @timeit
    @banner
    def Hook5(self, ice):
        "Prepare orientations for rigid water molecules."
        logger = getLogger()
        self.output["reppositions"] = ice.reppositions
        self.output["rotmatrices"]  = ice.rotmatrices


    @timeit
    @banner
    def Hook6(self, ice):
        "Arrange atoms of water and replacements."
        logger = getLogger()
        self.output["atoms"]        = ice.atoms


    @timeit
    @banner
    def Hook7(self, ice):
        "Arrange all atoms."
        logger = getLogger()
        self.output["atoms"]        = ice.atoms
