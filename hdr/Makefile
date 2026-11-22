CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: hdr

hdr: src/hdr.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: hdr
	./hdr

clean:
	rm -f hdr
