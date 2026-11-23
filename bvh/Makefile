CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra -Werror -ffp-contract=off

all: bvh

bvh: src/bvh.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: bvh
	./bvh

render: bvh
	./bvh render bvh.ppm

clean:
	rm -f bvh bvh.ppm
