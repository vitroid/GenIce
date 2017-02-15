#!/usr/bin/env python

from CifFile import ReadCif #CifFile
import logging


cf = ReadCif('RHO.cif')
logging.getLogger().info(cf)
