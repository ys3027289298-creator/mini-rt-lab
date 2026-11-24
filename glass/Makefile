CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: moire selftest

moire: src/moire.c src/glass.c src/glass.h
	$(CC) $(CFLAGS) src/moire.c src/glass.c -o $@ -lm

selftest: src/selftest.c src/glass.c src/glass.h
	$(CC) $(CFLAGS) src/selftest.c src/glass.c -o $@ -lm

test: selftest moire
	./selftest
	mkdir -p out
	./moire

figures: test
	python3 -I tools/plot.py out

clean:
	rm -f moire selftest
	rm -rf out

.PHONY: all test clean figures
