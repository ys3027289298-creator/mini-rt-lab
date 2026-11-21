CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: bvh qbvh

bvh: src/bvh.c
	$(CC) $(CFLAGS) $< -o $@ -lm

qbvh: src/qbvh.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: bvh qbvh
	./bvh
	./qbvh

clean:
	rm -f bvh qbvh
