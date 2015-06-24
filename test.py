#!/usr/bin/env python3

from subprocess import check_output
import sys

ices = "0 i 1c 1h 3 4 5 6 7 12 16 CS1 TS1 HS1 C0-II sTprime".split()

mode = "1"
if len(sys.argv) > 1:
    mode = sys.argv[1]

if mode == "2":
    for ice in ices:
        print("Ice",ice)
        result = check_output(["./genice", "-r","1","1","1", ice]).decode('utf-8')
        lines = result.split("\n")
        print(lines[-10:-1])
    sys.exit(0)
    


for ice in ices:
    for f in "mgdeq":
        result = check_output(["./genice", "-f", f, ice]).decode('utf-8')
        lines = result.split("\n")
        print(lines[-10:-1])
sys.exit(0)
