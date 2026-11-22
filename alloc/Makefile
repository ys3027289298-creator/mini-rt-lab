CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: alloc

alloc: src/alloc.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: alloc
	./alloc

clean:
	rm -f alloc
