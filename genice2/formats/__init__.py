# coding: utf-8

import sys


class Format():
    """
    Base class for Format()
    """

    def __init__(self, **kwargs):
        self.output = ""

    def hooks(self):
        return {}

    def dump(self):
        """
        This method returns self.output.
        """
        return self.output
