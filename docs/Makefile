build: parse
	bundle exec jekyll build --verbose --incremental
serve:
	cd ..; bundle exec jekyll serve --incremental --watch
%.png: %.svg
	inkscape $< -o $@
%x222.svg: Makefile
	 ../genice.x zeolite[$*] -f svg[shadow:rotatex=5:rotatey=5] -r 2 2 2 > $@
