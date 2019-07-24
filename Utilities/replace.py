#!/usr/bin/env python

# utility command

from __future__ import print_function
import sys
import string
tag  = sys.argv[1]
for line in open(sys.argv[2]):
    loc = line.find(tag)
    if loc >= 0:
        for newline in sys.stdin:
            print(" "*loc + newline, end="")
    else:
        print(line, end="")

