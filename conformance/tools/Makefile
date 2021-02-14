all: rms flipG192 mld

rms: rms.c tinywavein_c.h
	gcc -Wall -O2 -lm $< -o $@

flipG192: flipG192.c
	gcc -Wall -O2 -lm $< -o $@

mld:
	make -C mld_src/build config=release
	cp mld_src/bin/mld .
