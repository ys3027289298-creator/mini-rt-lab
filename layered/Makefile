CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: layered

layered: src/layered.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: layered
	./layered

clean:
	rm -f layered
