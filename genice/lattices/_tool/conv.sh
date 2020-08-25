#!/bin/bash
for i in *.py
do
    if [ ! -L $i ]
    then
	fgrep cellvectors $i || ( python convert.py `basename $i .py` > @; mv $i $i.old; mv @ $i)
    fi
done
