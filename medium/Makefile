CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: medium

medium: src/medium.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: medium
	./medium

clean:
	rm -f medium
