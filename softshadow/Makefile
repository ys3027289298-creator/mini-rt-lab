CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: softshadow

softshadow: src/softshadow.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: softshadow
	./softshadow

clean:
	rm -f softshadow
