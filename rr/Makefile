CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: rr

rr: src/rr.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: rr
	./rr

clean:
	rm -f rr
