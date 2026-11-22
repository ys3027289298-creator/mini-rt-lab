CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: density

density: src/density.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: density
	./density

clean:
	rm -f density
