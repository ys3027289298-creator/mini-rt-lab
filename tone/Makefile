CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: tone

tone: src/tone.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: tone
	./tone

clean:
	rm -f tone

.PHONY: all test clean
