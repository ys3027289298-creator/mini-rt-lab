CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: bvh packet

bvh: src/bvh.c
	$(CC) $(CFLAGS) $< -o $@ -lm

packet: src/packet.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: bvh packet
	./bvh
	./packet

clean:
	rm -f bvh packet
