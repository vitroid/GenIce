# coding: utf-8
import numpy as np
from logging import getLogger
import genice.molecules

class Molecule(genice.molecules.Molecule):
    def __init__(self):
    #placeholder for empty cage

        self.sites  = []
        self.labels = []
        self.name   = "empty"
