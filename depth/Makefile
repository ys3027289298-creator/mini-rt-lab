CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: depth

depth: src/depth.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: depth
	./depth

clean:
	rm -f depth
