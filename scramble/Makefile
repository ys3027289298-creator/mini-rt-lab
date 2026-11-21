CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: scramble

scramble: src/scramble.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: scramble
	./scramble

bench: scramble
	./scramble bench

dump: scramble
	./scramble dump

clean:
	rm -f scramble err_reuse.pgm err_full.pgm

.PHONY: all test bench dump clean
