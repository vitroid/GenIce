LIBPATH=/usr/local/share/genice
BINPATH=/usr/local/bin

all:
	echo Hello.
prepare:
	pip3 install numpy networkx
install: genice-install
	install genice-install $(BINPATH)/genice
	install -d $(LIBPATH)
	install -d $(LIBPATH)/lib
	install lib/*.py $(LIBPATH)/lib
	install -d $(LIBPATH)/Lattices
	install Lattices/*.py $(LIBPATH)/Lattices
	install -d $(LIBPATH)/Molecules
	install Molecules/*.py $(LIBPATH)/Molecules
genice-install:	genice
	sed -e 's@%%LIBPATH%%@$(LIBPATH)@' $< > $@
%: temp_%
	./genice -h | python3 lib/replace.py %%usage%% "    " $< > $@
test:
	i=23; while `python3 ./genice --rep 2 2 2 1c -f d --debug -s $$i > @`; do echo $$i; i=`expr $$i + 1`;done
distclean:
	-rm *.scad *.yap @*
	-rm -rf build dist
	-rm -rf GenIce.egg-info

