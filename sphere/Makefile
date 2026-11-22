CC ?= gcc
# -ffp-contract=off: the naive solver must not silently gain FMA accuracy,
# the stable solver calls fmaf explicitly.
CFLAGS ?= -O2 -std=c11 -Wall -Wextra -ffp-contract=off

OUT ?= ../artifacts/sphere

all: sphere

sphere: src/sphere.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: sphere
	mkdir -p $(OUT)
	./sphere $(OUT)

clean:
	rm -f sphere
