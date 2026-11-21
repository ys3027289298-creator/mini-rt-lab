CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: volume

volume: src/volume.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: volume
	./volume

clean:
	rm -f volume
