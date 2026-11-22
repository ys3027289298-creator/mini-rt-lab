CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: tonemap

tonemap: src/tonemap.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: tonemap
	./tonemap

clean:
	rm -f tonemap
