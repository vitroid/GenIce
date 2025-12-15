# coding: utf-8
"""
Centers of mass of water molecule
"""

from logging import getLogger
from genice2.decorators import timeit, banner
import genice2.formats


class Format(genice2.formats.Format):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def hooks(self):
        return {7: self.Hook7}

    @timeit
    @banner
    def Hook7(self, ice):
        "Do nothing."
        self.output = ""
