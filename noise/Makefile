CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: noise

noise: src/noise.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: noise
	./noise

clean:
	rm -f noise
