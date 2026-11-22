CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: uvdist

uvdist: src/uvdist.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: uvdist
	./uvdist --selftest

export: uvdist
	rm -rf out && mkdir -p out
	./uvdist --export out

clean:
	rm -f uvdist
	rm -rf out

.PHONY: all test export clean
