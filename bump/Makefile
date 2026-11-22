CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: bump

bump: src/bump.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: bump
	./bump

clean:
	rm -f bump
