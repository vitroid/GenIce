# coding: utf-8
import numpy as np
from logging import getLogger
import genice.molecules.one

class Molecule(genice.molecules.one.Molecule):
    def __init__(self):
        super().__init__(label="G15", name="C15")
