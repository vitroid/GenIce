# coding: utf-8
import genice3.molecule.one

desc = {"usage": "No options available.", "brief": "A united-atom methane model."}


class Molecule(genice3.molecule.one.Molecule):
    def __init__(self):
        super().__init__(label="Me", name="Me")
