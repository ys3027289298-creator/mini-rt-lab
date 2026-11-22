CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: converge

converge: src/converge.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: converge
	./converge --selftest

report: converge
	./converge

clean:
	rm -f converge
