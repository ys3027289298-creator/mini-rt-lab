CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra
RNG_CFLAGS := $(CFLAGS) -DRT_UNIFIED_RNG -I../shared

all: bvh

# Unified-source build (default): all randomness flows through shared/rng.c
bvh: src/bvh.c ../shared/rng.c ../shared/rng.h
	$(CC) $(RNG_CFLAGS) src/bvh.c ../shared/rng.c -o $@ -lm

# Legacy build: no RT_UNIFIED_RNG, original code paths, byte-identical output
bvh-legacy: src/bvh.c
	$(CC) $(CFLAGS) src/bvh.c -o $@ -lm

verify_image: verify_image.c src/bvh.c ../shared/rng.c
	$(CC) $(RNG_CFLAGS) -Dmain=bvh_selftest_main -Wno-unused-function verify_image.c ../shared/rng.c -o $@ -lm

test: bvh
	./bvh

test-legacy: bvh-legacy
	./bvh-legacy

clean:
	rm -f bvh bvh-legacy verify_image
