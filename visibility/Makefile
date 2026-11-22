CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: visibility

visibility: src/visibility.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: visibility
	./visibility

clean:
	rm -f visibility
