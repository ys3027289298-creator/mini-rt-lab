CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: onb

onb: src/onb.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: onb
	./onb selftest

sweep: onb
	mkdir -p out
	./onb sweep out/sweep.csv
	./onb flip out/flip.csv

plot: sweep
	python3 tools/plot_onb.py out/sweep.csv out/flip.csv out/onb_curves.png

clean:
	rm -rf onb out
