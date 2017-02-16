all:
	echo Hello.
%: temp_%
	./genice -h | python3 lib/replace.py %%usage%% "    " $< > $@
test:
	i=23; while `python3 ./genice --rep 2 2 2 1c -f d --debug -s $$i > @`; do echo $$i; i=`expr $$i + 1`;done
distclean:
	-rm *.scad *.yap @*
	-rm -rf build dist
	-rm -rf GenIce.egg-info

