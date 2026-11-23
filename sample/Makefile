CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra
RNG_CFLAGS := $(CFLAGS) -DRT_UNIFIED_RNG -I../shared

all: sample

# Unified-source build (default): all randomness flows through shared/rng.c
sample: src/sample.c ../shared/rng.c ../shared/rng.h
	$(CC) $(RNG_CFLAGS) src/sample.c ../shared/rng.c -o $@ -lm

# Legacy build: no RT_UNIFIED_RNG, original code paths, byte-identical output
sample-legacy: src/sample.c
	$(CC) $(CFLAGS) src/sample.c -o $@ -lm

verify_hist: verify_hist.c src/sample.c ../shared/rng.c
	$(CC) $(RNG_CFLAGS) -Dmain=sample_selftest_main -Wno-unused-function verify_hist.c ../shared/rng.c -o $@ -lm

test: sample
	./sample

test-legacy: sample-legacy
	./sample-legacy

clean:
	rm -f sample sample-legacy verify_hist
