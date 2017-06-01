#!/bin/bash

work=sort_tmp
mkdir $work
awk '/^G/{$1=""}{print}' $1 | sort > $work/a
awk '/^G/{$1=""}{print}' $2 | sort > $work/b
diff $work/a $work/b
