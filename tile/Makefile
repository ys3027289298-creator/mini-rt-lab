CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: tile

tile: src/tile.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: tile
	./tile

clean:
	rm -f tile
