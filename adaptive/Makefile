CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: adaptive

adaptive: src/adaptive.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: adaptive
	./adaptive

clean:
	rm -f adaptive error_se.pgm samples_n.pgm
