CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: envmap

envmap: src/envmap.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: envmap
	./envmap

clean:
	rm -f envmap
