%: temp_%
	./genice -h | python3 lib/replace.py %%usage%% "    " $< > $@

