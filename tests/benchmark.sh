\rm bench.out
for R in 2 4 8 16 24 32
do
  \rm @$R
  for rep in 1 2 3 4 5
  do
    ../genice.x ice1c -r $R $R $R --debug 2>> @$R > /dev/null
  done
  awk '/main: .* ms$/{s+=$5;n++}END{print 8*R*R*R,s/n}' R=$R @$R >> bench.out
done
cat bench.out
