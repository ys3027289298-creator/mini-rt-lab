CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: stattest

stattest: src/stattest.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: stattest
	./stattest

clean:
	rm -f stattest
