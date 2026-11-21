CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra -Werror
LDLIBS ?= -lm

all: intersect

intersect: src/intersect.c
	$(CC) $(CFLAGS) $< -o $@ $(LDLIBS)

test: intersect
	./intersect

clean:
	rm -f intersect

.PHONY: all test clean
