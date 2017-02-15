#!/usr/bin/env python3

# utility command

import sys
import string
tmp  = sys.argv[1]
by = sys.stdin.readlines()
sep = sys.argv[2]
by = sep.join(by)
by = sep+by
orig = open(sys.argv[3]).readlines()
orig = "".join(orig)
orig = orig.replace(tmp, by)
print(orig, end="")

