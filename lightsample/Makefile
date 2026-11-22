CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: lightsample

lightsample: src/lightsample.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: lightsample
	./lightsample

clean:
	rm -f lightsample
