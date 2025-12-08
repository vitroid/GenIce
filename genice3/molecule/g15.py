import numpy as np
import genice3.molecule.one


class Molecule(genice3.molecule.one.Molecule):
    def __init__(self):
        super().__init__(label="G15", name="C15")

