CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: offset

offset: src/offset.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: offset
	./offset

clean:
	rm -f offset
