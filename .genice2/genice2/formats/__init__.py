# coding: utf-8

import sys


class Format:
    """
    Base class for Format()
    """

    def __init__(self, **kwargs):
        raise NotImplementedError

    def dump(self, genice, file):
        """
        This method returns self.output.
        """
        raise NotImplementedError
