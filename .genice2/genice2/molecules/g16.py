# coding: utf-8
import numpy as np
from logging import getLogger
import genice2.molecules.one


class Molecule(genice2.molecules.one.Molecule):
    def __init__(self):
        super().__init__(label="G16", name="C16")
