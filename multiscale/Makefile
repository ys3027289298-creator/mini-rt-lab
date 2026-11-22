CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: multiscale

multiscale: src/multiscale.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: multiscale
	./multiscale

clean:
	rm -f multiscale
