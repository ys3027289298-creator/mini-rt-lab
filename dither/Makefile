CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: dither

dither: src/dither.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: dither
	./dither

clean:
	rm -f dither
