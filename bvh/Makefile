CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: bvh

bvh: src/bvh.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: bvh
	./bvh

clean:
	rm -f bvh
