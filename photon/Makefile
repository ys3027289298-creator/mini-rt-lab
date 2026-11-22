CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: photon

photon: src/photon.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: photon
	./photon

clean:
	rm -f photon
