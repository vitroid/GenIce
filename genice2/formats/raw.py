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
            1:self.hook1,
            2:self.hook2,
            3:self.hook3,
            4:self.hook4,
            5:self.hook5,
            6:self.hook6,
            7:self.hook7,
            }
        active = dict()
        for i in self.stages:
            if i in hook_funcs:
                active[i] = hook_funcs[i]
        return active


    def hook1(self, ice):
        logger = getLogger()
        logger.info("Hook1: Replicate water molecules to make a repeated cell.")
        self.output["reppositions"] = ice.reppositions
        self.output["repcell"]      = ice.repcell.mat
        self.output["repcagetype"]  = ice.repcagetype
        self.output["repcagepos"]   = ice.repcagepos
        self.output["cagetypes"]    = ice.cagetypes
        logger.info("Hook1: end.")


    def hook2(self, ice):
        logger = getLogger()
        logger.info("Hook2: Make a random graph and replicate.")
        self.output["dopants"]      = ice.dopants
        self.output["groups"]       = ice.groups
        self.output["filled_cages"] = ice.filled_cages
        self.output["graph"]        = ice.graph
        logger.info("Hook2: end.")


    def hook3(self, ice):
        logger = getLogger()
        logger.info("Hook3: Make a true ice graph.")
        self.output["graph"]        = ice.graph
        logger.info("Hook3: end.")


    def hook4(self, ice):
        logger = getLogger()
        logger.info("Hook4: Depolarize.")
        self.output["spacegraph"]   = ice.spacegraph
        logger.info("Hook4: end.")


    def hook5(self, ice):
        logger = getLogger()
        logger.info("Hook5: Prepare orientations for rigid water molecules.")
        self.output["reppositions"] = ice.reppositions
        self.output["rotmatrices"]  = ice.rotmatrices
        logger.info("Hook5: end.")


    def hook6(self, ice):
        logger = getLogger()
        logger.info("Hook6: Arrange atoms of water and replacements.")
        self.output["atoms"]        = ice.atoms
        logger.info("Hook6: end.")


    def hook7(self, ice):
        logger = getLogger()
        logger.info("Hook7: Arrange all atoms.")
        self.output["atoms"]        = ice.atoms
        logger.info("Hook7: end.")
