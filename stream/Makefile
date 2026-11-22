CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: stream

stream: src/stream.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: stream
	./stream

clean:
	rm -f stream
