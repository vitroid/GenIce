#!/usr/local/Cellar/gnuplot/5.2.7_1/bin/gnuplot

set xlabel "Number of molecules"
set ylabel "Time / msec"
set log xy
pl"bench.out" w linesp, "../../GenIce/tests/bench.out" w linesp, x
#    EOF
