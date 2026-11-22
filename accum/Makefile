CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: accum

accum: src/accum.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: accum
	./accum

clean:
	rm -f accum
