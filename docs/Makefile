build: parse
	bundle exec jekyll build --verbose --incremental
serve:
	cd ..; bundle exec jekyll serve --incremental --watch
%.png: %.svg
	inkscape $< -o $@
%x222.svg: Makefile
	../genice.x zeolite[$*] -f svg[shadow:rotatex=5:rotatey=5] -r 2 2 2 > $@
SOD.svg: Makefile
	../genice.x SOD -r 2 2 2 -f svg[polygon:shadow:rotatex=5:rotatey=-100]  > $@
#	../genice.x SOD -r 2 2 2 -f svg[polygon:shadow:rotatex=35:rotatey=-45]  > $@
FAU.svg: Makefile
	../genice.x FAU -r 2 2 2 -f svg[polygon:shadow:rotatex=5:rotatey=-100]  > $@
#	../genice.x FAU -r 2 2 2 -f svg[polygon:shadow:rotatex=35:rotatey=-45]  > $@
2xFAU.svg: Makefile
	../genice.x xFAU[2] -r 2 2 2 -f svg[polygon:shadow:rotatex=5:rotatey=-100]  > $@
#	../genice.x xFAU[2] -f svg[polygon:shadow:rotatex=35:rotatey=-45]  > $@
4xFAU.svg: Makefile
	../genice.x xFAU[4] -r 2 2 2 -f svg[polygon:shadow:rotatex=5:rotatey=-100]  > $@
#	../genice.x xFAU[4] -f svg[polygon:shadow:rotatex=35:rotatey=-45]  > $@
