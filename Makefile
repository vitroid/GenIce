%: temp_%
	./genice -h | python3 lib/replace.py %%usage%% "    " $< > $@
test:
	i=23; while ` ./genice --rep 2 2 2 1c -f d --debug -s $$i > @`; do echo $$i; i=`expr $$i + 1`;done

