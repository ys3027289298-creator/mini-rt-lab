CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra
RNG_CFLAGS := $(CFLAGS) -DRT_UNIFIED_RNG -I../shared

all: material

# Unified-source build (default): all randomness flows through shared/rng.c
material: src/material.c ../shared/rng.c ../shared/rng.h
	$(CC) $(RNG_CFLAGS) src/material.c ../shared/rng.c -o $@ -lm

# Legacy build: no RT_UNIFIED_RNG, original code paths, byte-identical output
material-legacy: src/material.c
	$(CC) $(CFLAGS) src/material.c -o $@ -lm

verify_variance: verify_variance.c src/material.c ../shared/rng.c
	$(CC) $(RNG_CFLAGS) -Dmain=material_selftest_main -Wno-unused-function verify_variance.c ../shared/rng.c -o $@ -lm

test: material
	./material

test-legacy: material-legacy
	./material-legacy

clean:
	rm -f material material-legacy verify_variance
