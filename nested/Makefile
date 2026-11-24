CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: nested

nested: src/nested.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: nested
	./nested selftest

clean:
	rm -f nested renders/*.ppm
