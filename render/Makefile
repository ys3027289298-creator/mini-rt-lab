CC ?= gcc
# -ffp-contract=off keeps the compiler from fusing multiply+add into FMA
# (which would make a pixel differ between builds/hosts); no -ffast-math
# so reassociation and NaN/Inf assumptions stay off. Determinism depends on
# these flags, and on -pthread supplying a real barrier-capable libc.
CFLAGS ?= -O2 -std=c11 -Wall -Wextra -pthread -ffp-contract=off

all: render

render: src/render.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: render
	./render

clean:
	rm -f render

.PHONY: all test clean
