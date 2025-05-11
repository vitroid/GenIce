# 改造したあとのbranchで実行する。
test:
	make -f Makefile.rules test

 
Makefile.rules: $(wildcard ../../genice2/*.py) mkmk.py
	python mkmk.py > $@

# 比較元(autotest-ref)のbranchで実行してファイルを準備する。
prepare: Makefile.rules
	make -f Makefile.rules -k prepare 

clean:
	-rm *.diff

distclean: clean
	-rm @*
	-rm *._KG
	-rm *._ringstat
	-rm *.cif
	-rm *.cif2
	-rm *.com
	-rm *.digraph
	-rm *.euler
	-rm *.exmol
	-rm *.exyz
	-rm *.graph
	-rm *.gromacs
	-rm *.mdview
	-rm *.output
	-rm *.povray
	-rm *.python
	-rm *.quaternion
	-rm *.rcom
	-rm *.rings
	-rm *.towhee
	-rm *.xyz
	-rm *.yaplot
	-rm Makefile.rules
	-rm *.diff

distclean: clean
	-rm @*
	-rm *._KG
	-rm *._ringstat
	-rm *.cif
	-rm *.cif2
	-rm *.com
	-rm *.digraph
	-rm *.euler
	-rm *.exmol
	-rm *.exyz
	-rm *.graph
	-rm *.gromacs
	-rm *.mdview
	-rm *.output
	-rm *.povray
	-rm *.python
	-rm *.quaternion
	-rm *.rcom
	-rm *.rings
	-rm *.towhee
	-rm *.xyz
	-rm *.yaplot
	-rm Makefile.rules
