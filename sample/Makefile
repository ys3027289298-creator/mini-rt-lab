CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: sample alloc_test

sample: src/sample.c
	$(CC) $(CFLAGS) $< -o $@ -lm

alloc_test: tests/alloc_test.c src/alloc.c src/alloc.h
	$(CC) $(CFLAGS) -Isrc tests/alloc_test.c src/alloc.c -o $@ -lm

test: sample alloc_test
	./sample
	./alloc_test

clean:
	rm -f sample alloc_test
