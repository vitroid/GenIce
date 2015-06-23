#!/bin/bash



for ice in 0 i 1c 1h 3 4 5 6 7 12 16 CS1 TS1 HS1 C0-II sTprime
do
    for f in m g d e q
    do
	./genice -f $f $ice | tail
    done
done
