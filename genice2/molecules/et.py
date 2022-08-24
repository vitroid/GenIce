# coding: utf-8
from logging import getLogger

import genice2.molecules.one

desc = {
    "usage": "No options available.",
    "brief": "A united-atom ethane model."
}


class Molecule(genice2.molecules.one.Molecule):
    def __init__(self):
        super().__init__(label="Et", name="Et")
