CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: penumbra

penumbra: src/penumbra.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: penumbra
	./penumbra

clean:
	rm -f penumbra
