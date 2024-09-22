# coding: utf-8
import numpy as np
from logging import getLogger
import genice2.molecules


class Molecule(genice2.molecules.Molecule):
    def __init__(self):
        # placeholder for empty cage

        self.sites_ = []
        self.labels_ = []
        self.name_ = "empty"
