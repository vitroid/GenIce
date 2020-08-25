# coding: utf-8

import sys

class Format():
    """
    Base class for Format()
    """
    def __init__(self, **kwargs):
        assert len(kwargs) == 0

    def hooks():
        return {}

    def dump(filename="", file=sys.stdout):
        pass # reserved for future
